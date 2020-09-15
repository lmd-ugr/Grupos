# Teoría-de-Grupos
Repositorio para TFG

Implementación de una librería en Python sobre teoría de Grupos de la asignatura ÁlgebraII del [Grado en Matemáticas](http://grados.ugr.es/matematicas/) de la [Universidad de Granada](http://www.ugr.es).

Extensión de la librería de José Luis Bueso y [pedritomelenas](https://github.com/pedritomelenas), basada en la librería de
[Naftali Harris](http://www.naftaliharris.com) y disponible en [Algebra-II](https://github.com/pedritomelenas/Algebra-II) y
[absalg](https://github.com/naftaliharris/Abstract-Algebra), respectivamente.

# Breve Descripción
- Un grupo es un conjunto G no vacío junto a una operación binaria * que verifica las tres propiedades de: asociatividad, existencia de elemento neutro y elemento inverso. Por ello, se usarán 3 ficheros.py principales:

- Set: Clase donde se añadirán todas las operaciones a nivel de conjunto.
- Function: Simula la operación binaria definida en nuestro grupo
- Group: Clase principal donde se realizarán las operaciones más importante:

- Adicionales: Se han añadido dos clases para representar las permutaciones y el grupo de los cuaternios
# Modificaciones

## Set.py
- Calcular el orden del conjunto
- Comprobar que el conjunto es finito
- Calcular los subconjuntos de un conjunto y los subconjuntos de tamaño n. Estas funciones son útiles de cara a la clase Grupo para así calcular los subgrupos y subgrupos de tamaño n con mayor facilidad.
- Además,se añaden las siguientes funciones que realizan operaciones a nivel de conjunto:
Unión, intersección, diferencia y diferencia simétrica.

## Function.py
- Método __str__: Lo implemento de nuevo para realizar una representación 
"bonita" y clara de nuestra función.
- El resto por ahora OK.

## Group.py

- Modifico Constructor:
	-Identidad: Puede ser que se le pase el elemento id erróneo por parámetro por lo que
	obligo a comprobar cual es el elemento neutro. 
	-Grupo abeliano. Borro la variable de clase y defino una función que compruebe
	si el grupo es abeliano (en vez de hacerlo en el constructor __init__)
	-Elimino la gran mayoría de las variables pasadas como parámetro en el constructor
	__init__ (check_associativity, check_inverses, identity) ya que el grupo debe
	cumplir estar 3 propiedades si o si, luego SIEMPRE SE OBLIGA a comprobar que el conjunto
	con la operación binaria cumplen las 3

- __str__: muestro además los elementos del grupo.


### class GroupElem:

- Inverse: La función que calcula la inversa de un elemento estaba definida en la clase Grupo en vez
de en la clase del elemento. Lo modifico.


- Siguiendo el libro ComputationalGroup Theory (libro verde), modifico 
__pow__ para realizar la exponenciación en O(log(n)) y también la
función order, que calcula el orden de los elementos del grupo con 
eficiencia O(log(n)^3) como mucho.


### class Group:

- Añado método Cardinality() (creo que ya estaba con el nombre de order() )


- Table -> Cayley_table  //pip install beautifultable
Elimino la implementación ENTERA ya que al usar HTML se utilizan librerías
que necesitan interfaz por lo que en la terminal no funcionaba el método. 
Además, el código era muy abstracto.
Simplicidad y optimización de código: En menos de la mitad de líneas realizo la misma funcionalidad
y se utiliza únicamente python. Funciona perfectamente en la terminal.
- De cara a funcionalidades próximas quizás es conveniente añadir representación
de la tabla con "letras".

- Añado función gens_cyclic_group: Obtiene los generadores del grupo.
Inicialmente estaba orientada para grupos cíclicos (phi Euler) pero
la he modificado para que funcione con cualquier grupo. 
Funciona pero aún no está terminada en su optimidad.



### Los siguientes métodos de subgrupos y normalidad los he modificado en su totalidad:

- El método is_subgroup(other) H <= G está mal, simplemente comprueba que
sea un subconjunto y que a*b (en G) = a*b (en H) \forall a,b \in G, es decir, que la operación
binaria se restringe. Pero NO comprueba la definición de subgrupo.
La he modificado y realizo lo siguiente:
- 1-Que sea un subconjunto H<=G (como a nivel de conjunto se ha implementado una función que
calcula los subconjuntos pues basta con comprobar que H es uno de ellos)
- 2- Si a,b \in H => a*b \in H 
- 3- Si a,b \in H => a*b^{-1} \in H


- Método all_subgroups(order): Recorro el contenedor de subconjuntos y devuelvo una 
lista con aquellos que sean subgrupos (llamo a la función is_subgroup(other) anterior.

- Método is_normalSubgroup(other): En caso de que sea un subgrupo compruebo 
que g*h*g**-1 pertenezca en el conjunto \forall g \in G , h in H


- Método all_normalSubgroups(order): Recorro el contenedor de subgrupos y a cada uno
de ellos le aplico is_normalSubgroup(other) devolviendo así todos los que
sean normales.

- Las funciones que devuelven todos los subgrupos o subgrupos normales llevan
además un parámetro "order" que nos dan la posibilidad de devolver aquellos subgrupos
de tamaño n=order indicado

- En todos estos métodos se devuelve una lista, quizás sea conveniente devolver un 
diccionario...? Por ahora no, funciona perfect.



## Permutation.py

- En una primera versión, las funciones que crean Simmetric y Alternative group las
había hecho en un método de clase. Sin embargo he preferido dejar el formato original y 
mantener esas funciones en el archivo Group.py

- Mágicamente la nueva tabla de Cayley funciona para todos los grupos :))

- He modificado en su totalidad la función __mul__ ya que estaba hecha bastante rara.
Modifico también la función __call__ para que funcione bien.

- Añado las funciones even_permutation y odd_permutation para saber si una
permutación es par o impar.

- Compruebo las nuevas funciones realizadas de la clase group (is_normal_group,
all_subgroups, is_normalSubgroup y all_normalSubgroups) con instancias
del grupo simétrico y alternado y funcionan perfect.

- Convendría quitar "=" en el __str__ y __repr__


## Quaternion.py

- Añado una Clase para encapsular los elementos de los cuaternios.

- Todos los operadores están incluídos así como funciones que nos permiten calcular
el conjugado, la norma, inverso, traza.

- Se puede instanciar un cuaternio a partir de una tupla o una letra (solo si es única):
Ejemplo: Quaternion(0,i,0,0) == Quaternion(letter="i")

- Gracias a lo anterior, ahora no es necesario indicar cómo se comporta el producto de los cuaternios pues todas 
estas operaciones ya están realizadas en la clase (__mul__)

- Se ha añadido un representación bonita de estos cuaternios (métodos __repr__ y __str__).
Se pueden realizar operaciones como i*k, k*i, j*k ...etc , la tabla de Cayley mostrada es correcta. 

-   >>> print(i*i == j*j == k*k == i*j*k == -1)
            'True' 


## Diedral.py

- Se ha creado una nueva clase para representar el grupo Diédrico de orden 2n. Funciona bien, mostrando todas las matrices
de rotación y reflexión de cada grupo.

- He tenido problemas a la hora de representar las rotaciones y reflexiones con matrices ya que no son hashables
y no se podía aplicar la operación binaria correctamente. He tenido que usar tuplas que al fin y al cabo representan
lo mismo.

- La tabla de Cayley para las matrices/tuplas es algo feota por lo que añado otra representación mejorada con R (rotaciones)
y S (simetrías/reflexiones). Se usa un diccionario para esta representación. RO, R1,...RN, S0, S1,... SN

- Todas las representaciones son equivalentes. Para probarlo se puede aplicar la función 'is_isomorphic', que devuelve True.

- Pulse [aquí](https://github.com/lmd-ugr/Grupos/blob/master/test/test_dihedral.png) para ver un ejemplo de las llamadas a la tabla de
Cayley con las tres diferentes representaciones.


## ToddCoxeter.py 

- Algoritmo de Todd Coxeter para resolver el problema de la palabra mediante la enumeración de clases.

- Se puede elegir el grupo del directorio 'Groups', devolviendo una tabla de Cosets de G/H y su cardinal, que 
como bien sabemos, coincide con el índice de G:H 

- Realizo otro pequeño algoritmo para obtener los generadores del grupo y, a partir de esos, se le da estructura de grupo.

- Funciona correctamente. Se puede dar un grupo por presentación y el mismo grupo dando los elementos. Al aplicar
is_isomorphic devuelve True.

## NEXT ToDo

- Decidirse por una función generate y generators buena.
- Modificar __pow__ de todas las clases? puede ser conveniente usar la función 
del libro ComputationalGroup theory por la eficiencia.
- Posible adición posterior de una representación matricial de los cuaternios.
- Repasar all_subgroups y all_normalSubgroups()
- Usar una clase para el grupo Diédrico (¡hecho!) y demás grupos que se vayan programando.
- Añadir operaciones al grupo diédrico (?)
- Quizás es conveniendo implementar un order() en cada grupo.

