# Mortalidad materna en Córdoba, Colombia, 2015 a 2024

Datos agregados y código para reproducir los resultados del estudio observacional de base poblacional sobre la mortalidad materna en el departamento de Córdoba entre 2015 y 2024. Este repositorio acompaña al manuscrito del mismo título.

## Disponibilidad de datos

Las fuentes primarias son públicas y oficiales: las estadísticas vitales del DANE (defunciones no fetales y nacimientos del registro civil). Este repositorio **no redistribuye datos a nivel de registro individual**. Para evitar la reidentificación de casos raros, en particular las muertes maternas en municipios pequeños, y en cumplimiento de la Ley 1581 de 2012 de protección de datos, solo se publican tablas agregadas y derivadas. La base a nivel de registro de las 162 muertes maternas tempranas puede solicitarse al autor de correspondencia, con fines de verificación, bajo un acuerdo de uso de datos.

## Estructura

- `Codigo/` — scripts de análisis y de figuras.
- `Datos derivados/` — bases agregadas (CSV) y la geometría municipal (GeoJSON).
- `Figuras/` — figuras del cuerpo del artículo y del material suplementario (PDF y PNG).
- `Salidas/` — resultados que produce el pipeline; se generan al ejecutarlo.

El diccionario de las columnas de cada base está en `DICCIONARIO_variables.md`.

## Cómo reproducir

1. Instalar las dependencias (Python 3.10 o superior):

   ```
   pip install -r requirements.txt
   ```

2. Ejecutar el pipeline desde la raíz del repositorio:

   ```
   python "Codigo/pipeline.py"
   ```

   Escribe `Salidas/resultados.txt` y un CSV por bloque de análisis. Los bloques que dependen de la base a nivel de registro (perfil sociodemográfico, ocurrencia por subregión) se omiten con un aviso cuando esa base no está presente.

3. Reproducir las figuras:

   ```
   python "Codigo/Figuras del articulo.py"
   ```

4. Comprobaciones de consistencia del paquete:

   ```
   python "Codigo/Autoauditoria estadistica.py"
   ```

5. Verificación independiente en R (magnitud, tendencia y análisis espacial):

   ```
   Rscript "Codigo/Analisis en R.R"
   ```

El archivo `Codigo/Codigo reproducible mortalidad.py` reconstruye las bases derivadas a partir de los microdatos originales del DANE. Esos microdatos no se incluyen aquí (ver *Disponibilidad de datos*) y no son necesarios para verificar los análisis agregados.

## Semillas

Permutaciones espaciales: 42. Simulación de potencia: 7, reiniciada en cada escenario, de modo que cada valor de potencia es reproducible por sí solo.
