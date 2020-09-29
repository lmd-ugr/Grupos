# Teoría-de-Grupos
Repositorio para TFG

Implementación de una librería en Python sobre teoría de Grupos de la asignatura Álgebra II del [Grado en Matemáticas](http://grados.ugr.es/matematicas/) de la [Universidad de Granada](http://www.ugr.es).

Extensión de la librería de José Luis Bueso y [pedritomelenas](https://github.com/pedritomelenas), basada en la librería de
[Naftali Harris](http://www.naftaliharris.com) y disponible en [Algebra-II](https://github.com/pedritomelenas/Algebra-II) y
[absalg](https://github.com/naftaliharris/Abstract-Algebra), respectivamente.

# Breve Descripción
- Un grupo es un conjunto G no vacío junto a una operación binaria * que verifica los siguientes axiomas: asociatividad, existencia de elemento neutro y existencia de elemento inverso para cada elemento. Por ello, se usarán 3 ficheros.py principales:

- Set: Clase donde se añadirán todas las operaciones a nivel de conjunto.
- Function: Simula la operación binaria definida en nuestro grupo
- Group: Clase principal donde se realizarán las operaciones más importantes:

- Ficheros Adicionales: Se han añadido dos clases para representar las permutaciones y el grupo de los cuaternios
# Modificaciones

## Set.py
Se han añadido bastantes métodos para que se puedan realizar operaciones a nivel de conjunto:
- Unión, diferencia, intersección, producto cartesiano y diferencia simétrica.

- `cardinality`, `is_finite`: Métodos para calcular la cardinalidad del conjunto y comprobar si este es finito.

- `subsets`, `subsets_n`: Métodos para calcular los subconjuntos de un conjunto y subconjunto de tamaño *n* . Estas dos métodos serán de mucha utilidad en la clase Grupo ya que simplificará mucho las operaciones.

## Function.py
- Se ha mantenido en tu totalidad el formato original, a excepción del operador `__str__` que muestra ahora la función de manera clara y precisa.


## Permutation.py

- Es la clase que construye el grupo simétrico y alternado. Para una mejor comprensión, creo el archivo *Permutation.py*  y añado toda la implementación existente.

- `__mul__`: Se modifica y se simplifica el producto de dos permutaciones.

- `__call__`: Este método fallaba cuando se llamaba con la imagen de *n* (longitud de la permutación). Se soluciona este error.

- `__even_permutation__`, `__odd_permutation__`: Métodos para calcular si una permutación es par o impar.


## ToddCoxeter.py

El algoritmo de **Todd Coxeter** es un algoritmo que resuelve el *problema de palabras* (*word problem*) para un grupo **G** mediante la enumeración de clases del grupo cociente **G/H**, donde **H** es un subgrupo de **G**.

La descripción del algoritmo se puede encontrar en la memoria del proyecto, la implementación en *ToddCoxeter.py* y un tutorial de su uso en *linkJupyter*. 

- `readGroup`: Implementación de una función que nos ayudará a leer los grupos por ficheros. Por orden, se leeran los generadores del grupo **G**, sus relaciones y los generadores del subgrupo **H**. En el directorio *Group* se proporcionaran ejemplos de grupos estudiados.

- `CosetEnumeration`: Método principal para llamar al algoritmo y obtener la tabla de clases de **G/H**.

- `schreier_graph`: Método que calcula el grafo de schreier resultante a partir de la tabla de clases laterales obtenidas del método anterior `CosetEnumeration`.

- `getGenerators`: A partir del grafo de Schreier no es difícil calcular el número de elementos del grupo. Para ello, calculamos sus generadores de forma recursiva usando este método.


Realizando operaciones con los generadores, se obtendrá el conjunto total de elementos al que le proporcionaremos estructura de grupo. 

- Se usará el *Teorema de Cayley* para representar cada grupo como grupo de permutaciones (usando las funcionalidades de *Permutation.py* y así usar una representación alternativa.



## Complex.py
Se ha realizado una implementación del grupo de las raíces n-ésimas de la unidad. Para ello, se ha 
implementado la clase número complejo junto a todos sus operadores que nos permiten sumar, restar, dividir, multiplicar...etc.

- `print_roots(roots)`: Función que representa las raíces *roots* pasadas como parámetro en el plano complejo.



## Quaternion.py

Se realiza una implementación de los números cuaternios en el archivo *Quaternion.py*.

- Se han implementado los principales operadores para trabajar con ellos, es decir, representarlos, sumar, restar, multiplicar, dividir, entre otros: `__repr__`, `__str__`, `__call__`, `__add__`, `__iadd__`, `__sub__`, `__eq__`, `__mull__`, `__div__` ...etc. Cabe destacar la importancia del operador `__mul__` ya que gracias a esta sobrecarga ya no hace falta indicar la tabla de multiplicar a la hora de crear el grupo.

- `conjugate`: Calcular el conjugado de un cuaternio.

- `norm`: Norma de un número cuaternio.

- `inverse`: Inverso.

- `trace`: Traza.

- Para crear un objeto se puede realizar de dos maneras, mediante su representación vectorial o indicando la letra en cuestión:
```python
>>> i = Quaternion(0,1,0,0)
>>> i2 = Quaternion(letter='i')
>>> i == i2
True
```

- Se puede probar de manera sencilla que las partes imaginarias verifican: 
```python
>>> i = Quaternion(0,1,0,0)
>>> j = Quaternion(0,0,1,0)
>>> j = Quaternion(0,0,0,1)
>>> i*i == j*j == k*k == i*j*k == -1
True
```

- Se verifican el resto de propiedades, como por ejemplo:
```python
>>> q = Quaternion(50,12,3,-9)
>>> r = Quaternion(-8,-2,2,32)
>>> (q*r).conjugate() == r.conjugate()*q.conjugate()
True
>>> (q*r).trace() == (r*q).trace()
True
```




## Diedral.py
Se ha realizado la implementación del grupo diédrico en el archivo *Dihedral.py*. Ahora, un grupo de orden *2n* estará formado por *n* simetrías y *n* rotaciones.

- Se han incorporado diferentes representaciones del grupo: *RS*, *tuple*, *permutations*
 -  RS: R0, R1,..., RN,  S0, S1,..., SN.
 - Tuple: representación mediante tuplas. (representan la matriz del movimiento)
 - Permutations: representación mediante permutaciones de *permutations.py*. 

Todas estas representaciones son equivalentes:
```python
>>> Dr = DihedralGroup(2, rep="RS")
>>> Dt = DihedralGroup(2, rep="tuple")
>>> Dp = DihedralGroup(2, rep="permutations")
>>> Dr.is_isomorphic(Dt)
True
>>> Dt.is_isomorphic(Dp)
True
>>> Dp.is_isomorphic(Dr)
True
```

- Todas las representaciones son equivalentes. Para probarlo se puede aplicar la función 'is_isomorphic', que devuelve True.

- En [aquí](https://github.com/lmd-ugr/Grupos/blob/master/test/test_dihedral.png) se puede ver un ejemplo de las llamadas a la tabla de Cayley con las tres diferentes representaciones.






## Group.py
En los archivos anteriores se han descrito las estructuras de diferentes grupos. Sin embargo, en *Group.py* se encuentran todos los métodos de nuestra librería. Se dividir a su vez en tres clases:
- Group:
- GroupElem:
- GroupAction:
- GroupHomomorphism:


### Class GroupElem
Representa un elemento de un grupo. Se ha añadido el método `inverse` a la implementación original. Este método calcula el inverso del elemento.

- Se ha añadido una implementación alternativa del operador `__pow__` para realizar la exponenciación en *O(log(n))*. 

- De igual modo, se ha añadido una implementación alternativa del método `order`, que calcula el orden de un elemento en *O(log(n)^3)* (como mucho).


### Class Group
- `is_abelian`: En la primera versión, se comprobaba si el grupo era abeliano en el constructor. Elimino la variable de clase y realizo esta comprobación en un método.


- `__str__` y `__repr__`: Se modifican para además mostrar los elementos del grupo (siempre que no teng un orden grande).

- `cosets`: Se optimiza y se simplifica el método.


- Table -> Cayley_table  //pip install beautifultable




## NEXT ToDo

- Repasar all_subgroups y all_normalSubgroups()
- Añadir operaciones al grupo diédrico (?)
- Arreglar la doble llamada a Todd Coxeter para obtener |G| y |H|.
- Añadir las instrucciones de uso al tutorial de Pedro.

