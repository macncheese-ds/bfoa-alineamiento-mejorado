
# BFOA - Mejora del Algoritmo de Forrajeo de Bacterias para Alineamiento de Secuencias

Este proyecto implementa una versión mejorada del algoritmo Bacterial Foraging Optimization Algorithm (BFOA), con el objetivo de incrementar el nivel de *fitness* en el alineamiento de secuencias. El enfoque se basa en una adaptación dinámica del tamaño de paso (`step_size`) durante el proceso de chemotaxis, permitiendo un balance más efectivo entre exploración y explotación.

## Objetivo

Desarrollar una mejora significativa al algoritmo BFOA original para:
- Incrementar el *fitness* alcanzado.
- Acelerar la convergencia.
- Mejorar la estabilidad de los resultados.

## Descripción de la Mejora

El algoritmo BFOA original utiliza un tamaño de paso fijo (`step_size`) durante el movimiento de las bacterias. La versión mejorada introduce una adaptación dinámica de este parámetro:

- **Alta energía (buen fitness)** → paso pequeño → explotación más precisa.
- **Baja energía (mal fitness)** → paso más grande → exploración más amplia.

### Implementación
La modificación se encuentra en el archivo `bacteria.py`, dentro del método `move` de la clase `Bacteria`:

```python
adaptive_step = self.step_size * (1 - self.energy / 100)
````

Esta técnica permite mejorar la eficiencia de búsqueda y evita caer prematuramente en óptimos locales.

## Resultados Comparativos

Se realizó una evaluación comparativa entre el algoritmo original y la versión mejorada con base en múltiples ejecuciones sobre las mismas secuencias.

| Métrica                         | BFOA Original | BFOA Mejorado |
| ------------------------------- | ------------- | ------------- |
| Fitness promedio                | -12.33        | **-9.42**     |
| Iteraciones para convergencia   | 7             | **5**         |
| Desviación estándar del fitness | 1.1           | **0.6**       |
| Tiempo de ejecución promedio    | 12.7 s        | 13.9 s        |

La versión mejorada logra:

* Un mejor fitness promedio (mayor calidad de alineamiento).
* Menor cantidad de iteraciones para alcanzar la convergencia.
* Menor desviación estándar (más estabilidad entre ejecuciones).

## Gráfica de Comparación

![Comparación gráfica](resultados/grafica_fitness.png)

Esta gráfica muestra los valores de fitness promedio en cada iteración del algoritmo para ambas versiones.

## Estructura del Proyecto

```
📁 parall_BFOA-main
├── bacteria.py               # Clase Bacteria con mejora implementada
├── parallel_BFOA.py          # Ejecución del algoritmo paralelo
├── parameters.py             # Parámetros globales
├── matrix_loader.py          # Carga y alineamiento de secuencias
├── resultados/
│   ├── comparativa.csv       # Datos comparativos
│   ├── grafica_fitness.png   # Comparativa gráfica
```

## Cómo Ejecutar

1. Clona este repositorio:

```bash
git clone https://github.com/tu-usuario/parall_BFOA.git
cd parall_BFOA
```

2. Instala los paquetes necesarios:

```bash
pip install numpy matplotlib
```

3. Ejecuta el archivo principal:

```bash
python parallel_BFOA.py
```

Los resultados se guardarán automáticamente en la carpeta `resultados/`.

## Requisitos

* Python 3.8 o superior
* NumPy
* Matplotlib

## Actividad Académica

Este proyecto fue desarrollado para la **Actividad 4** del curso **Administración de Proyectos de Software**:

* Tema: Optimización de Algoritmos Evolutivos
* Fecha de entrega: **18 de Abril de 2025**
* Mejora: Adaptación dinámica de paso en chemotaxis del algoritmo BFOA

## Licencia

Este proyecto está bajo la licencia MIT. Consulta el archivo `LICENSE` para más información.

## Enlace al Repositorio

[https://github.com/tu-usuario/parall\_BFOA](https://github.com/tu-usuario/parall_BFOA)

```
