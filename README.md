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

- Ficheros Adicionales: Se han añadido clases para representar el grupo de Permutaciones, Cuaternios, grupo de las raíces n-ésimas de la unidad y grupo Diédrico.

# Modificaciones

## Set.py
Se han añadido bastantes métodos para que se puedan realizar operaciones a nivel de conjunto:
- Unión, diferencia, intersección, producto cartesiano y diferencia simétrica.

- `cardinality`, `is_finite`: Métodos para calcular la cardinalidad del conjunto y comprobar si este es finito.

- `subsets`: Métodos para calcular los subconjuntos de un conjunto y subconjunto de tamaño *n* . Estas dos métodos serán de mucha utilidad en la clase Grupo ya que simplificará mucho las operaciones.

## Function.py
- Se ha mantenido en tu totalidad el formato original, a excepción del operador `__str__` que muestra ahora la función de manera clara y precisa.


```python
>>> S = Set({0,1,2})
>>> F = Function(S*S, S,lambda x: (x[0]+x[1])%3)
>>> print(F)
f((0, 1))=1   f((1, 2))=0   f((2, 1))=0
f((0, 0))=0   f((1, 1))=2   f((2, 0))=2
f((0, 2))=2   f((2, 2))=1   f((1, 0))=1
```


## Group.py
En los archivos anteriores se han descrito las estructuras de diferentes grupos. Sin embargo, en *Group.py* se encuentran todos los métodos de nuestra librería. Se divide a su vez en estas clases:
- Group:
- GroupElem:
- GroupAction:
- GroupHomomorphism:


### Class GroupElem
Representa un elemento de un grupo. Se ha añadido el método `inverse` a la implementación original. Este método calcula el inverso del elemento.

- Se ha añadido una implementación alternativa del operador `__pow__` para realizar la exponenciación en *O(log(n))*. 

- De igual modo, se ha añadido una implementación alternativa del método `order`, que calcula el orden de un elemento en *O(log(n)^3)* (como mucho).


### Class Group

- `__str__` y `__repr__`: Se modifican para además mostrar los elementos del grupo (siempre que no tengan un orden grande).

- `__init__`: Se modifica el constructor para poder definir grupos de varias formas:
 - Definición axiomatizada: Se comprueba que el par (Set,Function) pasado por argumento satisface los axiomas de grupo (asociatividad, identidad e inversos).

```python
>>> S = Set({0,1,2,3})
>>> F = Function(S*S, S,lambda x: (x[0]+x[1])%4)
>>> Z4 = Group(S,F)
>>> print(Z4)
Group with 4 elements: {0, 1, 2, 3}
```

 - Definición en términos de generadores y relatores. Sea un grupo $`\langle X \mid R \rangle`$. Se pasa por argumento el conjunto de generadores *X* y relaciones *R* que definen al grupo. El constructor se encarga de aplicar el **Algoritmo Todd Coxeter** y darle estructura de grupo de Permutaciones al grupo *G*.

```python
>>> gens = ['a']
>>> rels = ['aaaa'] #a^4=1
>>> G  = Group(gensG=gens, relsG=rels)
>>> print(G)
Group with 4 elements: {(), (1, 2, 3, 4), (1, 4, 3, 2), (1,3)(2, 4)}
```

Naturalmente, y aunque la forma de definir ambos grupos anteriores es distinta, son isomorfos:
```python
>>> G.is_isomorphic(Z4)
True
```

Por último, se añadirá una tercera forma de definir un grupo. SeaYunconjunto de elementos, entonces el grupo *G* se definirá como el grupo generado por$`\langle Y \rangle`$. En el siguiente ejemplo tomaremos un conjunto conuna única permutación, sin embargo, no exigimos que los elementossean permutaciones.
```python
>>> p = permutation((1,2,3,4))
>>> G = Group(elems=[p])
>>> print(G)
Group with 4 elements: {(), (1, 2, 3, 4), (1, 4, 3, 2), (1,3)(2, 4)}
```

- `is_abelian`: En la primera versión, se comprobaba si el grupo era abeliano en el constructor. Elimino la variable de clase y realizo esta comprobación en un método.

- `identity`: Del  mismo  modo  que  en *is_abelian*,  se  añade  un  nuevo  método para calcular la identidad del grupo.


- `cosets`:  Método  que  calcula  las  clases  laterales  de  un  grupo *G* sobre  un subgrupo *H*. Se optimiza y se simplifica.






## Permutation.py

- Es la clase que construye el grupo simétrico y alternado. Para una mejor comprensión, creo el archivo *Permutation.py*  y añado toda la implementación existente.

- `__mul__`: Se modifica y se simplifica el producto de dos permutaciones.

- `__call__`: Este método fallaba cuando se llamaba con la imagen de *n* (longitud de la permutación). Se soluciona este error.

- `__even_permutation__`, `__odd_permutation__`: Métodos para calcular si una permutación es par o impar.




## Complex.py
Se ha realizado una implementación del grupo de las raíces n-ésimas de la unidad. Para ello, se ha 
implementado la clase número complejo junto a todos sus operadores que nos permiten sumar, restar, dividir, multiplicar...etc.

- `plot(roots)`: Función que representa las raíces *roots* pasadas como parámetro en el plano complejo.



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


La función que se encarga de crear el grupo de los cuaternios es *QuaternionGroup, donde únicamente se le ha de pasar por argumento una de las dos representaciones siguientes:

```python
>>> Q = QuaternionGroup(rep="ijk")
>>> print(Q)
Group with 8 elements: { 1,  i,  j,  k,  -k,  -j,  -i,  -1}
```

```python
>>> Q2 = QuaternionGroup(rep="permutations")
>>> print(Q2)
Group with 8 elements: {(1, 4, 3, 2)(5, 7, 8, 6), (1, 7, 3,6)(2, 8, 4, 5), (1, 6, 3, 7)(2, 5, 4, 8), (1, 8, 3, 5)(2, 6,4, 7), (1, 2, 3, 4)(5, 6, 8, 7), (1, 5, 3, 8)(2, 7, 4, 6),(), (1, 3)(2, 4)(5, 8)(6, 7)}
```

```python
>>> Q.is_isomorphic(Q2)
True
```

- QuaternionGroupGeneralised(n): define el grupo generalizado de los cuaternios, con presentación:
$`Q_n = \langle a,b \mid a^n = b^2, a^{2n}=1,
b^{-1}ab=a^{-1} `$.



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




## ToddCoxeter.py

El  de **Algoritmo Todd Coxeter** es un algoritmo que resuelve el *Problema de Palabras* (*Word Problem*) para un grupo *G* mediante la enumeración de clases del grupo cociente *G/H* (a derechas), donde *H* es un subgrupo de *G*.

La descripción del algoritmo se puede encontrar en la memoria del proyecto, la implementación en [ToddCoxeter.py](https://github.com/lmd-ugr/Grupos/blob/master/ToddCoxeter.py) y un tutorial de su uso en [Jupyter](https://github.com/lmd-ugr/Grupos/blob/master/Tutorial.ipynb). 

- `readGroup`: Implementación de una función que nos ayudará a leer los grupos por ficheros. Por orden, se leeran los generadores del grupo *G*, sus relaciones y los generadores del subgrupo *H*. En el directorio *Group* se proporcionaran ejemplos de grupos estudiados.

- `CosetEnumeration`: Método principal para llamar al algoritmo y obtener la tabla de clases de *G/H*.

- `schreier_graph`: Método que calcula el grafo de schreier resultante a partir de la tabla de clases laterales obtenidas del método anterior `CosetEnumeration`.

- `getGenerators`: A partir del grafo de Schreier no es difícil calcular el número de elementos del grupo. Para ello, calculamos sus generadores de forma recursiva usando este método.


Realizando operaciones con los generadores, se obtendrá el conjunto total de elementos al que le proporcionaremos estructura de grupo. 

- Se usará el *Teorema de Cayley* para representar cada grupo como grupo de permutaciones (usando las funcionalidades de *Permutation.py* y así usar una representación alternativa.










