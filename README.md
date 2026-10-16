# Repositorio del Trabajo 1 de Series de Tiempo — ST8ZPT

**Repositorio:** https://github.com/acadavids-dot/st-2026-2-trabajo1-ST8ZPT
**Código del equipo:** ST8ZPT

## 1. Integrantes y aportes

| Integrante | Aportes al proyecto |
|---|---|
| Andrés Felipe Cadavid Sanguino | Responsable principal del desarrollo técnico del proyecto. Realizó gran parte de la implementación del código en R, el procesamiento y análisis de las series de tiempo, la evaluación y comparación de los métodos de pronóstico, y la elaboración del informe con sus resultados e interpretaciones. |
| Samuel Cárdenas García | Participó en tareas puntuales de apoyo y revisión del proyecto. Elaboró el archivo `README.md`, preparó la nota en PDF solicitada y diseñó las diapositivas de presentación. También apoyó la organización y presentación de los entregables del trabajo. |

## 2. Contenido del repositorio

| Archivo o carpeta | Contenido | Autor o responsable principal |
|---|---|---|
| `R/` | Scripts de R utilizados para implementar los procedimientos del análisis y los métodos de pronóstico. | Andrés Felipe Cadavid Sanguino |
| `analisis/correr-todo.R` | Script principal para ejecutar de forma integrada el análisis y generar los resultados del proyecto. | Andrés Felipe Cadavid Sanguino |
| `analisis/` | Archivos relacionados con la ejecución y organización de los análisis de las series. | Andrés Felipe Cadavid Sanguino |
| `datos/` | Datos utilizados en el análisis de las tres series de tiempo. | Andrés Felipe Cadavid Sanguino |
| `figs/` | Figuras y gráficos generados para el análisis y la presentación de los resultados. | Andrés Felipe Cadavid Sanguino, con apoyo puntual de Samuel Cárdenas García |
| `resultados/serie1.rds` | Resultados guardados correspondientes a la primera serie de tiempo. | Andrés Felipe Cadavid Sanguino |
| `resultados/serie2.rds` | Resultados guardados correspondientes a la segunda serie de tiempo. | Andrés Felipe Cadavid Sanguino |
| `resultados/serie3.rds` | Resultados guardados correspondientes a la tercera serie de tiempo. | Andrés Felipe Cadavid Sanguino |
| `informe/informe.qmd` | Archivo fuente del informe, elaborado en Quarto. | Andrés Felipe Cadavid Sanguino |
| `informe/informe.html` | Versión renderizada en HTML del informe del proyecto. | Andrés Felipe Cadavid Sanguino |
| `informe/secciones/` | Archivos Quarto que contienen las secciones del informe. | Andrés Felipe Cadavid Sanguino |
| `README.md` | Descripción del repositorio, instrucciones de ejecución, síntesis de resultados y declaración de uso de inteligencia artificial. | Samuel Cárdenas García |
| Nota en PDF solicitada | Documento complementario entregado como parte de los productos del trabajo. | Samuel Cárdenas García |
| Diapositivas de presentación | Material de apoyo para exponer los procedimientos y resultados del proyecto. | Samuel Cárdenas García |
| `.gitignore` | Archivos y carpetas que Git debe ignorar al registrar cambios. | Configuración del repositorio |
| `sesion-info.txt` | Archivo con información de la sesión de trabajo. | Según el contenido del archivo |

## 3. Instrucciones de ejecución

### Requisitos

- R instalado.
- Los paquetes de R requeridos por los scripts del proyecto.
- Quarto instalado para generar el informe en HTML.

### Ejecutar el análisis

Desde la carpeta raíz del repositorio, abrir R o RStudio y ejecutar:

```r
source("analisis/correr-todo.R")
```

También es posible medir el tiempo de ejecución con:

```r
system.time(source("analisis/correr-todo.R"))
```

En la ejecución realizada, se obtuvo el siguiente tiempo:

| Medida | Tiempo |
|---|---:|
| Tiempo de usuario | 122.85 segundos |
| Tiempo de sistema | 4.30 segundos |
| Tiempo transcurrido (`elapsed`) | 180.21 segundos |

El tiempo transcurrido fue de aproximadamente **3 minutos**. Este valor corresponde a la ejecución registrada y puede variar según el equipo y el entorno de R.

### Generar el informe en HTML

Desde una terminal ubicada en la carpeta raíz del repositorio, ejecutar:

```bash
quarto render informe/informe.qmd
```

Este comando renderiza el archivo fuente de Quarto y genera el informe en HTML. El tiempo de renderizado debe registrarse por separado si se requiere reportarlo.

## 4. Resumen de resultados

La evaluación de los métodos se realizó mediante pronósticos con origen móvil y comparación del RMSE frente a un método de referencia. En la tabla, cada RMSE se presenta como **RMSE del ganador / RMSE del referente**, seguido de la razón entre ambos. Una razón inferior a 1 indica que el método ganador obtuvo menor RMSE que el referente.

| Serie | Transformación seleccionada | Ganador en \(h=1\) | Parámetros en \(h=1\) | RMSE ganador / referente en \(h=1\) | Ganador en \(H\) | Parámetros en \(H\) | RMSE ganador / referente en \(H\) |
|---|---|---|---|---|---|---|---|
| 1. Llegadas de turistas | \(\nabla(1-B^{12})\log(Y_t)\) | Holt–Winters multiplicativo (`hw_mult`) | \(s=12,\ \alpha=0.5,\ \beta=0.1,\ \gamma=0.1\) | 2.96 / 6.77; razón = 0.437 | Holt–Winters multiplicativo (`hw_mult`) | \(s=12,\ \alpha=0.5,\ \beta=0.1,\ \gamma=0.1\) | 4.19 / 7.76; razón = 0.541 |
| 2. Habitaciones de hotel ocupadas | \(\nabla Y_t=(1-B)Y_t\) | Suavizamiento exponencial simple (`ses`) | \(\alpha=0.72\) | 6.24 / 6.49; razón = 0.962 | Tendencia cuadrática (`tend_cuad`) | Coeficientes estimados mediante el ajuste cuadrático | 8.31 / 12.10; razón = 0.688 |
| 3. Estudiantes matriculados | \(Y_t-\widehat{\beta}_0-\widehat{\beta}_1t\) | Tendencia lineal (`tend_lineal`) | Coeficientes estimados por mínimos cuadrados | 10.80 / 14.20; razón = 0.762 | Doble media móvil (`dmm`) | \(k=12\) | 10.40 / 14.50; razón = 0.712 |

**Interpretación general.** En las tres series, los métodos seleccionados obtuvieron un RMSE inferior al de sus respectivos referentes tanto para \(h=1\) como para el horizonte máximo \(H\). En la primera serie, Holt–Winters multiplicativo capturó la estacionalidad anual. En la segunda, el suavizamiento exponencial simple obtuvo el mejor resultado a un paso, mientras que la tendencia cuadrática se destacó en el horizonte máximo. En la tercera, la tendencia lineal fue la mejor opción a un paso y la doble media móvil obtuvo el menor RMSE para \(H=6\). Estas diferencias deben interpretarse teniendo en cuenta la variabilidad de los errores y la cercanía entre algunos métodos.

## 5. Declaración de uso de inteligencia artificial

Durante la elaboración del proyecto se utilizaron herramientas de inteligencia artificial como apoyo complementario al trabajo realizado por los integrantes. Su uso incluyó la depuración y revisión de código, la identificación de posibles errores de ejecución, el apoyo en ajustes y modificaciones de figuras, y la mejora de la redacción, organización y claridad de algunos apartados del informe y de la documentación del repositorio.

La inteligencia artificial se utilizó como herramienta de asistencia y no como sustituto de la revisión del equipo. Los resultados, procedimientos estadísticos y contenidos entregados deben contrastarse con las salidas obtenidas en R y con los criterios establecidos para el trabajo. La responsabilidad sobre el contenido final corresponde a los integrantes del equipo.