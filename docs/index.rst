vmath
-----

.. py:class:: Vector
.. py:class:: Quaternion
.. py:class:: Matrix

.. py:method:: vmath.vec(x: float, y: float, z: float) -> Vector
.. py:method:: vmath.quat(x: float, y: float, z: float, w: float) -> Quaternion
.. py:method:: vmath.mat(position: Vector, rotation: Quaternion, scale: Vector) -> Matrix

.. py:method:: vmath.rotate_x(angle: float) -> Quaternion
.. py:method:: vmath.rotate_y(angle: float) -> Quaternion
.. py:method:: vmath.rotate_z(angle: float) -> Quaternion

.. py:method:: vmath.rotate(axis: Vector, angle: float) -> Quaternion
.. py:method:: vmath.slerp(a: Quaternion, b: Quaternion, t: float) -> Quaternion

.. py:method:: vmath.random_axis() -> Vector
.. py:method:: vmath.random_rotation() -> Quaternion

Vector
------

.. py:method:: Vector.normal() -> Vector
.. py:method:: Vector.length() -> float
.. py:method:: Vector.pack() -> bytes

Quaternion
----------

.. py:method:: Quaternion.axis() -> Vector
.. py:method:: Quaternion.angle() -> float
.. py:method:: Quaternion.inverse() -> Quaternion
.. py:method:: Quaternion.pack() -> bytes

Matrix
------

.. py:method:: Matrix.inverse() -> Matrix
.. py:method:: Matrix.position() -> Vector
.. py:method:: Matrix.rotation() -> Quaternion
.. py:method:: Matrix.scale() -> Vector
.. py:method:: Matrix.pack() -> bytes

Scaling a Vector
----------------

.. py:method:: Vector * float -> Vector
    :noindex:

.. py:method:: float * Vector -> Vector
    :noindex:

.. py:method:: Vector * Vector -> Vector
    :noindex:

Rotate a Vector
---------------

.. py:method:: Quaternion * Vector -> Vector
    :noindex:

Rotate a Quaternion
-------------------

.. py:method:: Quaternion * Quaternion -> Quaternion
    :noindex:

Apply Transform
---------------

.. py:method:: Matrix * Matrix -> Matrix
    :noindex:

.. py:method:: Matrix * Vector -> Vector
    :noindex:
