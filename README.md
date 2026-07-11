# recortar_hojas

Recorta una ortofoto grande en un GeoTIFF por cada hoja de un índice `.ASC`.

Cada hoja del índice es un cuadrilátero girado (no un rectángulo). El programa escribe un
GeoTIFF por hoja, con el norte arriba y **sin remuestrear**: los píxeles salen tal cual están en
la ortofoto de origen. Lo que cae fuera del polígono de la hoja se graba transparente.

La ortofoto se lee por bloques de filas, así que **da igual lo grande que sea**: el programa gasta
siempre la misma memoria. Una ortofoto de 3 GB no da ningún problema.

---

## 1. Instalar Python

En Windows 10 y 11, desde la línea de órdenes:

```
winget install --id Python.Python.3.12 --scope machine
```

`--scope machine` lo instala para todos los usuarios del ordenador. Si prefiere instalarlo sólo
para usted, quite esa parte.

> Si `winget` no funciona en su ordenador, descargue el instalador de
> <https://www.python.org/downloads/windows/> y ejecútelo marcando la casilla
> **«Add python.exe to PATH»**.

**Cierre la ventana de la consola y abra otra**, para que Windows se entere de que Python está
instalado. Después compruebe que responde:

```
python --version
```

Debe contestar `Python 3.12.x` o superior.

---

## 2. Instalar las dependencias

El programa necesita **rasterio**, que trae consigo `numpy` y la biblioteca GDAL. Ocupa unos
134 MB en total.

```
python -m pip install rasterio
```

Compruebe que ha quedado bien instalado:

```
python -c "import rasterio; print(rasterio.__version__)"
```

Debe contestar con un número de versión (por ejemplo `1.5.0`) y nada más. Si dice
`ModuleNotFoundError`, la instalación no ha funcionado.

---

## 3. Uso

Forma general:

```
python recortar_hojas.py ORTOFOTO HOJAS EPSG [opciones]
```

Los tres primeros valores son **obligatorios** y van siempre en ese orden:

| Valor      | Qué es                                                                    |
| ---------- | ------------------------------------------------------------------------- |
| `ORTOFOTO` | La ortofoto de origen, en TIFF.                                            |
| `HOJAS`    | El índice de hojas, en `.ASC`.                                             |
| `EPSG`     | El código del sistema de referencia que se graba en la cabecera de cada hoja. |

Si el nombre de un fichero lleva espacios, hay que escribirlo **entre comillas**.

### Lo primero: mirar sin escribir nada

Antes de ponerse a recortar conviene ver qué hojas hay y lo que va a ocupar cada una. La opción
`--simular` enseña la lista y no toca el disco:

```
python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --simular
```

### Recortar de verdad

```
python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --salida hojas
```

Deja un fichero por hoja (`H1.tif`, `H2.tif`...) dentro de la carpeta `hojas`, que se crea sola si
no existe. Las hojas que ya estén hechas **se respetan**: para rehacerlas hay que añadir
`--sobrescribir`.

### Recortar sólo algunas hojas

```
python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --solo H2,H6 --sobrescribir
```

### Ver la ayuda

```
python recortar_hojas.py --ayuda
```

---

## 4. Opciones

| Opción                     | Para qué sirve                                                                     |
| -------------------------- | ---------------------------------------------------------------------------------- |
| `-s`, `--salida`           | Carpeta donde se dejan las hojas. Por omisión, `hojas`.                            |
| `--orden-ejes xy` \| `yx`  | Orden en que están grabadas las coordenadas de la ortofoto. Véase el apartado 5.    |
| `--blanco-transparente`    | Graba como transparente el blanco puro. **Activado por omisión.** Para desactivarlo, `--no-blanco-transparente`. |
| `--solo H2,H6`             | Recorta sólo las hojas que se indiquen, separadas por comas.                       |
| `--simular`                | Enseña la lista de hojas y su tamaño, sin escribir nada.                           |
| `--sobrescribir`           | Vuelve a generar las hojas cuyo fichero ya exista.                                 |
| `--omitir-vacias`          | No crea fichero para las hojas que caen fuera de la ortofoto.                      |
| `--plantilla-nombre`       | Nombre de cada fichero. Admite `{nombre}`, `{indice}` y `{orto}`. Por omisión, `{nombre}.tif`. |
| `-q`, `--silencioso`       | No enseña el porcentaje de avance.                                                 |

Hay además un grupo de **opciones avanzadas** (tamaño de tesela, memoria de GDAL...) que casi
nunca hace falta tocar. Salen todas con `--ayuda`.

### Que las hojas ocupen menos

Las hojas salen comprimidas en **DEFLATE**, que no pierde nada de calidad y lo lee cualquier
programa, por viejo que sea. Si le parecen muy grandes, hay dos maneras de reducirlas:

```
python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --compresion WEBP
```

**WEBP ocupa la mitad y no pierde absolutamente nada**: los píxeles salen idénticos, bit a bit.
La hoja H6, por ejemplo, pasa de 252 MB a 135 MB. Tarda unos tres minutos por hoja en vez de uno.
El único inconveniente es que WEBP dentro de un TIFF es un formato moderno: QGIS y ArcGIS
actuales lo leen sin problema, pero un programa antiguo puede no saber abrirlo. **Si va a
entregar las hojas a alguien, asegúrese antes de que su programa las abre.**

Si aún así necesita que ocupen mucho menos, se puede bajar la calidad de WEBP:

```
python recortar_hojas.py A-44.tif "Hojas 5km_A-44.asc" 25830 --compresion WEBP --calidad-webp 95
```

Con calidad 95 la hoja ocupa **siete veces menos** y la pérdida no se aprecia a simple vista. Pero
es una pérdida de verdad: **la calidad que se tira no se recupera nunca**. Piénselo dos veces si
las hojas son una entrega para un cliente.

> **JPEG no sirve.** Es la compresión que a todo el mundo se le ocurre primero, pero no sabe
> llevar la transparencia: la estropea. El programa se niega a usarla.

---

## 5. Dos cosas que conviene entender

### El código EPSG

El programa **no se fía** del sistema de referencia que traiga la ortofoto: graba en cada hoja el
que usted le indique en la línea de órdenes, y le dice por pantalla cuál traía el original.

Para la ortofoto de la A-44 el código correcto es **25830** (ETRS89 / UTM zona 30N). El fichero de
origen viene marcado como EPSG:3042, que es el mismo sistema pero con los ejes al revés, y hay
muchos programas que no lo entienden.

### El orden de las coordenadas (`--orden-ejes`)

Lo normal es que la ortofoto llegue **sin sistema de referencia**: en la cabecera sólo están la
coordenada de la primera esquina y el tamaño del píxel. Esas coordenadas van en orden **X,Y**
(este, norte), que es lo que el programa supone y **no hay que indicar nada**.

Algunos sistemas, como el EPSG:3042, declaran los ejes al revés (norte-este), y hay programas que
graban entonces la coordenada en orden **Y,X**. Sólo en ese caso hay que decírselo:

```
python recortar_hojas.py orto.tif hojas.asc 25830 --orden-ejes yx
```

El programa **no lo adivina**, pero sí avisa: si la ortofoto y las hojas no llegan a tocarse —que
es lo que pasa cuando el orden no es el bueno—, lo dice por pantalla y le sugiere probar el otro.
Si no avisa de nada, el orden es el correcto.

---

## 6. Si algo va mal

| Lo que dice la pantalla                              | Qué pasa                                                                 |
| ---------------------------------------------------- | ------------------------------------------------------------------------ |
| `ModuleNotFoundError: No module named 'rasterio'`    | Falta instalar rasterio. Véase el apartado 2.                            |
| `error: no se encuentra ...`                         | El nombre del fichero está mal escrito, o falta ponerlo entre comillas.  |
| `AVISO: ... la ortofoto y las hojas no coinciden`    | El orden de los ejes no es el bueno. Véase el apartado 5.                |
| Las hojas salen en blanco                            | Lo mismo: pruebe con `--orden-ejes yx`.                                  |

Las hojas ocupan bastante: cada una de las de 5 km de la A-44 pesa entre **260 y 410 MB**, y las
seis juntas pasan de **2 GB**. Asegúrese de que hay sitio en el disco antes de empezar.
# recortar_hojas
