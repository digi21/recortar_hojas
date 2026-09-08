#!/usr/bin/env python3
"""Recorta una ortofoto grande en un GeoTIFF por cada hoja descrita en un índice .ASC.

La ortofoto de origen se lee por bloques de filas, así que el consumo de memoria es
constante sea cual sea el tamaño del fichero. La imagen se relee una vez por hoja.

Formato .ASC
------------
Líneas (polígonos):
    C=<código> <n_puntos> <indicador>
    <x> <y> <z>            (repetido n_puntos veces)

Textos:
    T=<código> <n_puntos> <indicador>
    <altura_texto> <justificación> <rotación>
    <x> <y> <z>            (repetido n_puntos veces)
    <texto>

Cada hoja es un polígono más el texto colocado dentro de él, que da el nombre de la hoja.

Orden de los ejes
-----------------
Lo normal es que la ortofoto llegue sin sistema de referencia, o con uno desconocido: en
la cabecera del GeoTIFF sólo están la coordenada de la primera esquina y el tamaño del
píxel. Esas coordenadas van en orden X,Y (este, norte), que es lo que el programa supone.

Hay sistemas, como el EPSG:3042 («ETRS89 / UTM zone 30N (N-E)»), que declaran los ejes al
revés, norte-este, y algunos programas graban entonces la coordenada en orden Y,X. GDAL
nunca invierte nada al leer: se cree siempre que lo grabado es X,Y. Para esos ficheros
está la opción --orden-ejes yx, que hay que pedir a mano. El programa no lo adivina, pero
avisa si la ortofoto y el índice de hojas no llegan a tocarse, que es lo que pasa cuando
el orden elegido no es el bueno.

El blanco es el vacío
---------------------
La ortofoto original no tiene banda de transparencia: da por vacío el blanco puro
(255,255,255). Por eso, y salvo que se diga --no-blanco-transparente, todo píxel blanco
puro se graba transparente. Se ha comprobado en la ortofoto de trabajo que el blanco puro
sólo aparece en el fondo, nunca dentro de la imagen: ni las marcas viales de la autovía
llegan a saturar las tres bandas a 255.

Salida
------
Un GeoTIFF por hoja, norte arriba y sin remuestreo: la ventana es el rectángulo
envolvente del polígono de la hoja. Los píxeles que caen fuera del polígono (o fuera de
la extensión del origen) se escriben totalmente transparentes. El código EPSG indicado
se graba en la cabecera; el que trajera la ortofoto no se tiene en cuenta.

Con --relleno R G B, los píxeles transparentes que caen dentro del polígono se graban
con ese color y opacos; los de fuera del polígono siguen transparentes.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.enums import ColorInterp
from rasterio.features import geometry_mask
from rasterio.windows import Window
from rasterio.windows import transform as transformada_ventana

Punto = tuple[float, float]
Rectangulo = tuple[float, float, float, float]  # x_min, y_min, x_max, y_max


# --------------------------------------------------------------------------- #
# Traducción de los mensajes propios de argparse
# --------------------------------------------------------------------------- #

# argparse pasa sus textos por gettext, así que basta con sustituir la función que usa
# internamente para que toda su salida (ayuda y errores) salga en español.

_MENSAJES = {
    "usage: ": "uso: ",
    "positional arguments": "argumentos obligatorios",
    "options": "opciones",
    "show this help message and exit": "muestra esta ayuda y termina",
    "the following arguments are required: %s": "faltan estos argumentos obligatorios: %s",
    "unrecognized arguments: %s": "no se entienden estos argumentos: %s",
    "invalid %(type)s value: %(value)r": "el valor %(value)r no es un %(type)s válido",
    "argument %(argument_name)s: %(message)s": "en %(argument_name)s: %(message)s",
    "expected one argument": "se esperaba un valor",
    "expected at least one argument": "se esperaba al menos un valor",
    "expected %s argument": "se esperaba %s valor",
    "expected %s arguments": "se esperaban %s valores",
    "invalid choice: %(value)r (choose from %(choices)s)":
        "%(value)r no es una opción válida; elija entre %(choices)s",
    "ambiguous option: %(option)s could match %(matches)s":
        "la opción %(option)s es ambigua: puede ser %(matches)s",
    "not allowed with argument %s": "no se puede usar junto con %s",
    "ignored explicit argument %r": "se ignora el valor %r",
    "%(prog)s: error: %(message)s\n": "%(prog)s: error: %(message)s\n",
}


def _traducir(mensaje: str) -> str:
    return _MENSAJES.get(mensaje, mensaje)


def _traducir_plural(singular: str, plural: str, n: int) -> str:
    return _traducir(singular if n == 1 else plural)


argparse._ = _traducir            # type: ignore[attr-defined]
argparse.ngettext = _traducir_plural  # type: ignore[attr-defined]


def _entero(texto: str) -> int:
    return int(texto)


def _ruta(texto: str) -> Path:
    return Path(texto)


# argparse nombra el tipo con `type.__name__` al quejarse de un valor mal escrito.
_entero.__name__ = "número entero"
_ruta.__name__ = "nombre de fichero"


class _Ayuda(argparse.RawDescriptionHelpFormatter):
    """Muestra el valor por defecto de cada opción, en español, cuando es informativo."""

    def _get_help_string(self, accion: argparse.Action) -> str:
        ayuda = accion.help or ""
        # Se comparan identidades: `2 in (True,)` sería cierto, y ocultaría un valor real.
        sin_interes = any(
            accion.default is valor for valor in (None, False, True, argparse.SUPPRESS)
        )
        muestra_defecto = accion.option_strings and "%(default)" not in ayuda and not sin_interes
        return ayuda + " (por defecto: %(default)s)" if muestra_defecto else ayuda


# --------------------------------------------------------------------------- #
# Lectura del .ASC
# --------------------------------------------------------------------------- #

@dataclass(slots=True)
class Poligono:
    codigo: str
    anillo: list[Punto]  # cerrado, el primer punto se repite al final


@dataclass(slots=True)
class Texto:
    codigo: str
    punto: Punto
    valor: str


@dataclass(slots=True)
class Hoja:
    nombre: str
    anillo: list[Punto]


def _leer_lineas(ruta: Path) -> list[str]:
    """Lee el fichero probando UTF-8 primero y luego las codificaciones habituales."""
    datos = ruta.read_bytes()
    for codificacion in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return datos.decode(codificacion).splitlines()
        except UnicodeDecodeError:
            continue
    raise ValueError(f"no se puede descodificar {ruta}")


def _coordenadas(linea: str) -> Punto:
    partes = linea.split()
    if len(partes) < 2:
        raise ValueError(f"se esperaba una línea de coordenadas, se encontró {linea!r}")
    return float(partes[0]), float(partes[1])


def leer_asc(ruta: Path) -> tuple[list[Poligono], list[Texto]]:
    lineas = _leer_lineas(ruta)
    poligonos: list[Poligono] = []
    textos: list[Texto] = []
    i = 0

    while i < len(lineas):
        linea = lineas[i].strip()
        if not linea:
            i += 1
            continue

        cabecera = linea.split()
        if "=" not in cabecera[0]:
            raise ValueError(f"{ruta}:{i + 1}: se esperaba una cabecera de entidad, "
                             f"se encontró {linea!r}")

        tipo, _, codigo = cabecera[0].partition("=")
        tipo = tipo.upper()
        try:
            n_puntos = int(cabecera[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{ruta}:{i + 1}: número de puntos incorrecto "
                             f"en {linea!r}") from exc

        inicio = i
        i += 1
        if tipo == "T":
            i += 1  # altura, justificación, rotación

        if i + n_puntos > len(lineas):
            raise ValueError(f"{ruta}:{inicio + 1}: entidad cortada al final del fichero")
        puntos = [_coordenadas(lineas[i + k]) for k in range(n_puntos)]
        i += n_puntos

        if tipo == "C":
            poligonos.append(Poligono(codigo, _cerrar(puntos)))
        elif tipo == "T":
            valor = lineas[i].strip() if i < len(lineas) else ""
            i += 1
            textos.append(Texto(codigo, puntos[0], valor))
        else:
            print(f"aviso: {ruta}:{inicio + 1}: se ignora la entidad {tipo}=", file=sys.stderr)

    return poligonos, textos


def _cerrar(puntos: list[Punto]) -> list[Punto]:
    if len(puntos) < 3:
        raise ValueError(f"polígono con {len(puntos)} puntos")
    anillo = list(puntos)
    if anillo[0] != anillo[-1]:
        anillo.append(anillo[0])
    if len(anillo) < 4:
        raise ValueError("polígono degenerado")
    return anillo


def punto_en_anillo(punto: Punto, anillo: list[Punto]) -> bool:
    """Punto dentro de polígono, por lanzamiento de rayo sobre un anillo cerrado."""
    x, y = punto
    dentro = False
    for (x1, y1), (x2, y2) in zip(anillo, anillo[1:]):
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            dentro = not dentro
    return dentro


_ILEGALES = re.compile(r'[<>:"/\\|?*\x00-\x1f]')


def sanear(nombre: str) -> str:
    """Quita del nombre los caracteres que Windows no admite en un fichero."""
    limpio = _ILEGALES.sub("_", nombre).strip().rstrip(". ")
    return limpio or "sin_nombre"


def construir_hojas(poligonos: list[Poligono], textos: list[Texto]) -> list[Hoja]:
    """Empareja cada polígono con el texto colocado en su interior."""
    hojas: list[Hoja] = []
    usados: set[int] = set()

    for indice, poligono in enumerate(poligonos, start=1):
        candidatos = [
            j for j, texto in enumerate(textos)
            if j not in usados and punto_en_anillo(texto.punto, poligono.anillo)
        ]
        if candidatos:
            if len(candidatos) > 1:
                print(
                    f"aviso: hay {len(candidatos)} textos dentro del polígono n.º {indice}, "
                    f"se usa {textos[candidatos[0]].valor!r}",
                    file=sys.stderr,
                )
            elegido = candidatos[0]
            usados.add(elegido)
            nombre = textos[elegido].valor
        else:
            nombre = f"{poligono.codigo or 'hoja'}_{indice}"
            print(
                f"aviso: no hay ningún texto dentro del polígono n.º {indice}, "
                f"se le llama {nombre!r}",
                file=sys.stderr,
            )
        hojas.append(Hoja(sanear(nombre), poligono.anillo))

    return _desduplicar(hojas)


def _desduplicar(hojas: list[Hoja]) -> list[Hoja]:
    """Si dos hojas se llaman igual, añade un número al nombre de la segunda."""
    vistos: dict[str, int] = {}
    for hoja in hojas:
        cuenta = vistos.get(hoja.nombre, 0) + 1
        vistos[hoja.nombre] = cuenta
        if cuenta > 1:
            hoja.nombre = f"{hoja.nombre}_{cuenta}"
    return hojas


# --------------------------------------------------------------------------- #
# Orden de los ejes de la ortofoto
# --------------------------------------------------------------------------- #

def transformada_origen(origen: rasterio.DatasetReader, orden_ejes: str) -> Affine:
    """Georreferenciación de la ortofoto, con los ejes en el orden que se indique.

    GDAL da siempre por hecho que lo grabado en la cabecera es X,Y. Si el fichero lo trae
    al revés (norte, este), hay que deshacer el cambio: se intercambian la coordenada de
    la esquina y el tamaño de píxel de cada eje.
    """
    t = origen.transform
    if orden_ejes == "xy":
        return t
    if t.b or t.d:
        raise ValueError(
            "la ortofoto está girada (no tiene el norte arriba) y no se sabe cómo "
            "invertirle los ejes; hay que georreferenciarla de nuevo con otro programa"
        )
    # t.a es el ancho del píxel y t.e el alto, en negativo porque las filas van hacia el sur.
    return Affine(-t.e, t.b, t.f, t.d, -t.a, t.c)


_NOMBRE_SISTEMA = re.compile(r'(?:PROJCS|GEOGCS)\["([^"]+)"')


def sistema_declarado(origen: rasterio.DatasetReader) -> str:
    """Qué sistema de referencia dice traer la ortofoto. Sólo sirve para informar: el que
    se graba en las hojas es el que se pide por la línea de órdenes.
    """
    if not origen.crs:
        return "ninguno"
    wkt = origen.crs.to_wkt()
    nombre = coincidencia.group(1) if (coincidencia := _NOMBRE_SISTEMA.search(wkt)) else None
    # Con el umbral de confianza que trae de fábrica (70), to_epsg() no reconoce el
    # EPSG:3042 de la ortofoto de ejemplo, porque su datum viene en forma de conjunto.
    codigo = origen.crs.to_epsg(confidence_threshold=25)
    if codigo and nombre:
        return f"EPSG:{codigo} «{nombre}»"
    if codigo:
        return f"EPSG:{codigo}"
    return f"«{nombre}», sin código EPSG" if nombre else "desconocido"


def rectangulo_ortofoto(origen: rasterio.DatasetReader, transformada: Affine) -> Rectangulo:
    esquinas = [
        transformada * (columna, fila)
        for columna in (0, origen.width)
        for fila in (0, origen.height)
    ]
    xs, ys = zip(*esquinas)
    return min(xs), min(ys), max(xs), max(ys)


def rectangulo_hojas(hojas: list[Hoja]) -> Rectangulo:
    puntos = [punto for hoja in hojas for punto in hoja.anillo]
    xs, ys = zip(*puntos)
    return min(xs), min(ys), max(xs), max(ys)


def _solapan(uno: Rectangulo, otro: Rectangulo) -> bool:
    return (
        uno[0] < otro[2] and otro[0] < uno[2]  # en X
        and uno[1] < otro[3] and otro[1] < uno[3]  # en Y
    )


def avisar_si_no_solapan(
    origen: rasterio.DatasetReader,
    transformada: Affine,
    hojas: list[Hoja],
    orden_ejes: str,
) -> None:
    """Avisa si la ortofoto y el índice de hojas ni siquiera se tocan.

    Es lo que ocurre cuando el orden de los ejes elegido no es el bueno, y si no se dice
    el resultado son hojas enteras en blanco sin explicación. No se corrige nada por
    cuenta propia: sólo se avisa y se sugiere qué probar.
    """
    if _solapan(rectangulo_ortofoto(origen, transformada), rectangulo_hojas(hojas)):
        return

    otro = "yx" if orden_ejes == "xy" else "xy"
    print(
        f"AVISO: leyendo las coordenadas de la ortofoto en orden {orden_ejes.upper()}, la "
        f"ortofoto y las hojas no coinciden: no tienen ni un punto en común, así que "
        f"todas las hojas van a salir en blanco.\n"
        f"       Es lo que pasa cuando el orden de los ejes no es el bueno. "
        f"Pruebe a añadir «--orden-ejes {otro}».",
        file=sys.stderr,
    )


# --------------------------------------------------------------------------- #
# Recorte
# --------------------------------------------------------------------------- #

@dataclass(slots=True)
class Bandas:
    """Qué bandas del origen hay que leer y cómo se escriben."""
    color: list[int]        # índices (empezando en 1) de las bandas de color del origen
    alfa: int | None        # índice (empezando en 1) de la banda alfa del origen, si la hay

    @property
    def n_salida(self) -> int:
        return len(self.color) + 1


def distribucion_bandas(origen: rasterio.DatasetReader) -> Bandas:
    interpretaciones = list(origen.colorinterp)
    if ColorInterp.alpha in interpretaciones:
        alfa = interpretaciones.index(ColorInterp.alpha) + 1
        color = [i for i in range(1, origen.count + 1) if i != alfa]
    else:
        alfa = None
        color = list(range(1, origen.count + 1))
    if not color:
        raise ValueError("el origen no tiene bandas de color")
    return Bandas(color, alfa)


def valor_opaco(tipo_dato: str) -> int | float:
    """Valor de la banda alfa que significa «totalmente opaco» para ese tipo de dato."""
    return np.iinfo(tipo_dato).max if np.issubdtype(np.dtype(tipo_dato), np.integer) else 1


def ventana_hoja(transformada: Affine, anillo: list[Punto]) -> Window:
    """Rectángulo envolvente del anillo, en coordenadas de píxel del origen.

    La ventana puede salirse del origen; eso se resuelve al leer.
    """
    inversa = ~transformada
    columnas, filas = zip(*(inversa * punto for punto in anillo))
    col_ini = math.floor(min(columnas))
    fila_ini = math.floor(min(filas))
    ancho = math.ceil(max(columnas)) - col_ini
    alto = math.ceil(max(filas)) - fila_ini
    return Window(col_ini, fila_ini, ancho, alto)


def opciones_creacion(args: argparse.Namespace, bandas: Bandas) -> dict:
    compresion = args.compresion.upper()
    opciones = {
        "driver": "GTiff",
        "tiled": True,
        "blockxsize": args.tam_tesela,
        "blockysize": args.tam_tesela,
        "compress": compresion,
        "bigtiff": args.bigtiff,
        "sparse_ok": args.disperso,
        "num_threads": "ALL_CPUS",
        "photometric": "RGB" if len(bandas.color) == 3 else "MINISBLACK",
        "alpha": "YES",
    }
    if compresion in ("DEFLATE", "LZW", "ZSTD", "LERC_DEFLATE", "LERC_ZSTD"):
        opciones["predictor"] = args.predictor
    if compresion == "DEFLATE":
        opciones["zlevel"] = args.nivel_deflate
    if compresion == "WEBP":
        # La calidad 100 es la compresión sin pérdida de WEBP: los píxeles salen idénticos,
        # y aun así ocupa casi la mitad que DEFLATE. Por debajo de 100 hay pérdida.
        if args.calidad_webp >= 100:
            opciones["webp_lossless"] = True
        else:
            opciones["webp_lossless"] = False
            opciones["webp_level"] = args.calidad_webp
    return opciones


def recortar_hoja(
    origen: rasterio.DatasetReader,
    transformada: Affine,
    hoja: Hoja,
    bandas: Bandas,
    crs: CRS,
    destino: Path,
    args: argparse.Namespace,
) -> bool:
    """Escribe una hoja. Devuelve False si la hoja no solapa con el origen."""
    ventana = ventana_hoja(transformada, hoja.anillo)
    ancho_salida, alto_salida = int(ventana.width), int(ventana.height)
    transformada_salida = transformada_ventana(ventana, transformada)
    col_ini, fila_ini = int(ventana.col_off), int(ventana.row_off)

    solapa = (
        col_ini < origen.width and fila_ini < origen.height
        and col_ini + ancho_salida > 0 and fila_ini + alto_salida > 0
    )
    if not solapa and args.omitir_vacias:
        return False

    geometria = {"type": "Polygon", "coordinates": [hoja.anillo]}
    tipo_dato = origen.dtypes[0]
    opaco = valor_opaco(tipo_dato)
    indices_lectura = bandas.color + ([bandas.alfa] if bandas.alfa else [])
    relleno = None if args.relleno is None else np.array(args.relleno, dtype=tipo_dato)

    tesela = args.tam_tesela
    filas_bloque = max(tesela, math.ceil(args.filas_bloque / tesela) * tesela)

    perfil = opciones_creacion(args, bandas)
    perfil.update(
        width=ancho_salida, height=alto_salida, count=bandas.n_salida,
        dtype=tipo_dato, crs=crs, transform=transformada_salida,
    )

    with rasterio.open(destino, "w", **perfil) as destino_ds:
        destino_ds.colorinterp = _interpretacion_color(len(bandas.color))

        for fila_salida in range(0, alto_salida, filas_bloque):
            filas = min(filas_bloque, alto_salida - fila_salida)
            ventana_bloque = Window(0, fila_salida, ancho_salida, filas)
            dentro = geometry_mask(
                [geometria],
                out_shape=(filas, ancho_salida),
                transform=transformada_ventana(ventana_bloque, transformada_salida),
                all_touched=args.todo_tocado,
                invert=True,
            )

            columnas = np.flatnonzero(dentro.any(axis=0))
            if columnas.size == 0:
                _progreso(hoja.nombre, fila_salida + filas, alto_salida, args.silencioso)
                continue  # en estas filas no hay nada de la hoja: no se tocan las teselas

            # Se ajusta la franja escrita a la rejilla de teselas de salida para que GDAL
            # nunca tenga que releer una tesela escrita a medias.
            primera = (int(columnas[0]) // tesela) * tesela
            ultima = min(ancho_salida, math.ceil((int(columnas[-1]) + 1) / tesela) * tesela)
            dentro = dentro[:, primera:ultima]
            franja = ultima - primera

            buffer = np.zeros((bandas.n_salida, filas, franja), dtype=tipo_dato)

            fila_origen = fila_ini + fila_salida
            col_origen = col_ini + primera
            f0, f1 = max(fila_origen, 0), min(fila_origen + filas, origen.height)
            c0, c1 = max(col_origen, 0), min(col_origen + franja, origen.width)

            if f1 > f0 and c1 > c0:
                datos = origen.read(
                    indexes=indices_lectura,
                    window=Window(c0, f0, c1 - c0, f1 - f0),
                )
                df, dc = f0 - fila_origen, c0 - col_origen
                filas_leidas = slice(df, df + (f1 - f0))
                cols_leidas = slice(dc, dc + (c1 - c0))
                buffer[:len(indices_lectura), filas_leidas, cols_leidas] = datos
                if bandas.alfa is None:
                    # El origen no tiene alfa: todo lo que se ha leído es opaco.
                    buffer[-1, filas_leidas, cols_leidas] = opaco

                if args.blanco_transparente:
                    # En la ortofoto original el blanco es el color del vacío, y no hay
                    # banda de transparencia. Se pasa esa convención a la banda alfa.
                    blanco = (buffer[:len(bandas.color)] == opaco).all(axis=0)
                    buffer[-1][blanco] = 0

            # Fuera del polígono: totalmente transparente, y el color a cero para que
            # las teselas se compriman hasta casi nada.
            buffer *= dentro

            if relleno is not None:
                # Ya está a cero todo lo de fuera del polígono, así que un alfa a cero
                # dentro de `dentro` es un hueco de la hoja: el vacío de la ortofoto o
                # la parte de la hoja a la que la ortofoto no llega.
                hueco = dentro & (buffer[-1] == 0)
                buffer[:len(bandas.color), hueco] = relleno[:, None]
                buffer[-1][hueco] = opaco

            destino_ds.write(buffer, window=Window(primera, fila_salida, franja, filas))
            _progreso(hoja.nombre, fila_salida + filas, alto_salida, args.silencioso)

    if not args.silencioso:
        print(file=sys.stderr)
    return True


def _interpretacion_color(n_color: int) -> list[ColorInterp]:
    if n_color == 3:
        return [ColorInterp.red, ColorInterp.green, ColorInterp.blue, ColorInterp.alpha]
    if n_color == 1:
        return [ColorInterp.gray, ColorInterp.alpha]
    return [ColorInterp.undefined] * n_color + [ColorInterp.alpha]


def _progreso(nombre: str, hechas: int, total: int, silencioso: bool) -> None:
    if silencioso:
        return
    print(f"\r  {nombre}: {100.0 * hechas / total:5.1f}%", end="", file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# Punto de entrada
# --------------------------------------------------------------------------- #

_EJEMPLOS = """\
Ejemplos:

  Ver las hojas y el tamaño que tendrá cada una, sin escribir nada:
    python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --simular

  Recortar todas las hojas en la carpeta «hojas»:
    python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 -s hojas

  Recortar sólo dos hojas, pisando los ficheros que ya existan:
    python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --solo H2,H6 --sobrescribir

  Ortofoto cuyas coordenadas están grabadas en orden norte-este (Y,X):
    python recortar_hojas.py orto.tif hojas.asc 25830 --orden-ejes yx

  Rellenar de negro los huecos de dentro de cada hoja:
    python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --relleno 0 0 0
"""


def analizar_argumentos(argv: list[str] | None = None) -> argparse.Namespace:
    analizador = argparse.ArgumentParser(
        description="Recorta una ortofoto en un GeoTIFF por cada hoja de un índice .ASC.",
        epilog=_EJEMPLOS,
        formatter_class=_Ayuda,
        add_help=False,
    )
    analizador.add_argument("-a", "--ayuda", "-h", "--help", action="help",
                            help="muestra esta ayuda y termina")
    analizador.add_argument("ortofoto", type=_ruta, help="ortofoto de origen (TIFF)")
    analizador.add_argument("hojas", type=_ruta, help="índice de hojas (.ASC)")
    analizador.add_argument("epsg", type=_entero,
                            help="código EPSG que se graba en la cabecera del GeoTIFF")
    analizador.add_argument("-s", "--salida", type=_ruta, default=Path("hojas"),
                            help="carpeta donde se dejan las hojas recortadas")
    analizador.add_argument("--orden-ejes", choices=("xy", "yx"), default="xy",
                            help="orden en que están grabadas las coordenadas de la "
                                 "ortofoto: xy = este, norte (lo habitual); yx = norte, "
                                 "este (algunos ficheros en EPSG:3042 y demás sistemas "
                                 "de eje invertido)")
    analizador.add_argument("--plantilla-nombre", default="{nombre}.tif",
                            help="nombre de cada fichero; admite {nombre}, {indice} y {orto}")
    analizador.add_argument("--solo", default=None,
                            help="nombres de hoja, separados por comas, que hay que recortar; "
                                 "por omisión se recortan todas")
    analizador.add_argument("--blanco-transparente", action=argparse.BooleanOptionalAction,
                            default=True,
                            help="da por vacío el blanco puro (255,255,255) de la ortofoto y "
                                 "lo graba como transparente. Con --no-blanco-transparente el "
                                 "blanco se conserva tal cual")
    analizador.add_argument("--relleno", type=_entero, nargs=3, metavar=("R", "G", "B"),
                            default=None,
                            help="rellena con este color, y deja opacos, los píxeles "
                                 "transparentes que caigan dentro de la hoja: el vacío de la "
                                 "ortofoto y la parte de la hoja a la que la ortofoto no "
                                 "llega. Fuera de la hoja se sigue grabando transparente")
    analizador.add_argument("--simular", action="store_true",
                            help="enseña la lista de hojas y su tamaño, sin escribir nada")
    analizador.add_argument("--sobrescribir", action="store_true",
                            help="vuelve a generar las hojas cuyo fichero ya exista")
    analizador.add_argument("--omitir-vacias", action="store_true",
                            help="no crea fichero para las hojas que caen fuera de la ortofoto")
    analizador.add_argument("-q", "--silencioso", action="store_true",
                            help="no enseña el porcentaje de avance")

    avanzadas = analizador.add_argument_group("opciones avanzadas")
    avanzadas.add_argument("--compresion", default="DEFLATE",
                           help="compresión de la salida: DEFLATE, WEBP, ZSTD, LZW o NONE. "
                                "WEBP ocupa casi la mitad que DEFLATE sin perder un solo bit, "
                                "pero los programas antiguos no saben leerlo")
    avanzadas.add_argument("--calidad-webp", type=_entero, default=100,
                           help="calidad de WEBP, de 1 a 100. 100 es sin pérdida: los píxeles "
                                "salen idénticos. Por debajo de 100 la hoja ocupa mucho menos, "
                                "pero la imagen se degrada y eso ya no se deshace")
    avanzadas.add_argument("--nivel-deflate", type=_entero, default=6,
                           help="nivel de compresión de DEFLATE, de 1 a 9")
    avanzadas.add_argument("--predictor", type=_entero, default=2, choices=(1, 2, 3),
                           help="predictor del TIFF; 2 = diferencias horizontales")
    avanzadas.add_argument("--bigtiff", default="IF_SAFER",
                           choices=("YES", "NO", "IF_NEEDED", "IF_SAFER"),
                           help="cuándo usar el formato BigTIFF, necesario por encima de 4 GB")
    avanzadas.add_argument("--tam-tesela", type=_entero, default=512,
                           help="lado de la tesela de salida, en píxeles")
    avanzadas.add_argument("--filas-bloque", type=_entero, default=512,
                           help="filas que se leen y escriben de una vez; a más filas, "
                                "más memoria. Se redondea a una fila entera de teselas")
    avanzadas.add_argument("--cache-gdal", type=_entero, default=512,
                           help="memoria de trabajo de GDAL, en MB")
    avanzadas.add_argument("--todo-tocado", action=argparse.BooleanOptionalAction, default=True,
                           help="incluye los píxeles que sólo roza el borde del polígono")
    avanzadas.add_argument("--disperso", action="store_true",
                           help="deja sin escribir las teselas que caen fuera de la hoja; "
                                "ocupa algo menos, pero hay programas que no leen los "
                                "GeoTIFF dispersos")
    return analizador.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = analizar_argumentos(argv)

    if not args.ortofoto.is_file():
        print(f"error: no se encuentra {args.ortofoto}", file=sys.stderr)
        return 1
    if not args.hojas.is_file():
        print(f"error: no se encuentra {args.hojas}", file=sys.stderr)
        return 1

    try:
        crs = CRS.from_epsg(args.epsg)
    except Exception as exc:
        print(f"error: el código EPSG {args.epsg} no es válido: {exc}", file=sys.stderr)
        return 1

    if args.compresion.upper() == "JPEG":
        # Comprobado: JPEG no sabe llevar la banda de transparencia, la altera. Las hojas
        # saldrían con la transparencia estropeada, y sin avisar.
        print("error: JPEG no vale, porque estropea la transparencia de las hojas. "
              "Para que ocupen menos, use «--compresion WEBP».", file=sys.stderr)
        return 1
    if not 1 <= args.calidad_webp <= 100:
        print("error: --calidad-webp va de 1 a 100", file=sys.stderr)
        return 1

    poligonos, textos = leer_asc(args.hojas)
    hojas = construir_hojas(poligonos, textos)
    if not hojas:
        print(f"error: no hay ningún polígono en {args.hojas}", file=sys.stderr)
        return 1

    if args.solo:
        pedidas = {sanear(nombre.strip()) for nombre in args.solo.split(",")}
        desconocidas = pedidas - {hoja.nombre for hoja in hojas}
        if desconocidas:
            print(f"error: estas hojas no están en el índice: {', '.join(sorted(desconocidas))}",
                  file=sys.stderr)
            return 1
        hojas = [hoja for hoja in hojas if hoja.nombre in pedidas]

    entorno = rasterio.Env(GDAL_CACHEMAX=args.cache_gdal, GDAL_NUM_THREADS="ALL_CPUS")
    with entorno, rasterio.open(args.ortofoto) as origen:
        bandas = distribucion_bandas(origen)

        if args.relleno is not None:
            if len(bandas.color) != 3:
                cuantas = "1 banda" if len(bandas.color) == 1 else f"{len(bandas.color)} bandas"
                print(f"error: --relleno da un color R,G,B y esta ortofoto tiene {cuantas} "
                      f"de color, no 3", file=sys.stderr)
                return 1
            maximo = valor_opaco(origen.dtypes[0])
            if not all(0 <= componente <= maximo for componente in args.relleno):
                print(f"error: los valores de --relleno van de 0 a {maximo} en esta ortofoto",
                      file=sys.stderr)
                return 1

        try:
            transformada = transformada_origen(origen, args.orden_ejes)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1

        traia = sistema_declarado(origen)
        ejes = "este, norte" if args.orden_ejes == "xy" else "norte, este"
        print(
            f"{args.ortofoto.name}: {origen.width} x {origen.height} píxeles, "
            f"{origen.count} bandas ({len(bandas.color)} de color y "
            f"{'1 de transparencia' if bandas.alfa else 'ninguna de transparencia'}), "
            f"píxel de {abs(transformada.a)} m\n"
            f"  sistema que traía la ortofoto: {traia}; se graba EPSG:{args.epsg}\n"
            f"  coordenadas leídas en orden {args.orden_ejes.upper()} ({ejes})",
            file=sys.stderr,
        )
        avisar_si_no_solapan(origen, transformada, hojas, args.orden_ejes)

        if args.simular:
            for indice, hoja in enumerate(hojas, start=1):
                ventana = ventana_hoja(transformada, hoja.anillo)
                crudo = int(ventana.width) * int(ventana.height) * bandas.n_salida
                print(
                    f"  {indice:3d}  {hoja.nombre:<20s} "
                    f"{int(ventana.width):7d} x {int(ventana.height):7d} píxeles "
                    f"({crudo / 2**30:6.2f} GiB sin comprimir)"
                )
            return 0

        args.salida.mkdir(parents=True, exist_ok=True)
        escritas = omitidas = 0

        for indice, hoja in enumerate(hojas, start=1):
            destino = args.salida / args.plantilla_nombre.format(
                nombre=hoja.nombre, indice=indice, orto=args.ortofoto.stem
            )
            if destino.exists() and not args.sobrescribir:
                print(f"[{indice}/{len(hojas)}] {destino.name}: ya existe, se omite",
                      file=sys.stderr)
                omitidas += 1
                continue

            print(f"[{indice}/{len(hojas)}] {destino.name}", file=sys.stderr)
            if recortar_hoja(origen, transformada, hoja, bandas, crs, destino, args):
                escritas += 1
            else:
                print("  queda fuera de la ortofoto, se omite", file=sys.stderr)
                omitidas += 1

    print(f"terminado: {escritas} hojas escritas, {omitidas} omitidas, en {args.salida}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
