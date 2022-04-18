from typing import Iterable


class Vector:
    x: float
    y: float
    z: float
    def normal(self) -> Vector: ...
    def length(self) -> float: ...
    def pack(self) -> bytes: ...
    def __neg__(self) -> Vector: ...
    def __add__(self, other: Vector) -> Vector: ...
    def __sub__(self, other: Vector) -> Vector: ...
    def __mul__(self, other: Vector | float) -> Vector: ...
    def __rmul__(self, other: Matrix | float) -> Vector: ...
    def __iter__(self) -> Iterable[float]: ...
    def __getitem__(self, key: int) -> float: ...


class Quaternion:
    x: float
    y: float
    z: float
    w: float
    def axis(self) -> Vector: ...
    def angle(self) -> float: ...
    def inverse(self) -> Quaternion: ...
    def pack(self) -> bytes: ...
    def __mul__(self, other: Quaternion) -> Quaternion: ...
    def __iter__(self) -> Iterable[float]: ...
    def __getitem__(self, key: int) -> float: ...


class Matrix:
    def inverse(self) -> Matrix: ...
    def position(self) -> Vector: ...
    def rotation(self) -> Quaternion: ...
    def scale(self) -> Vector: ...
    def pack(self) -> bytes: ...
    def __mul__(self, other: Matrix) -> Matrix: ...
    def __iter__(self) -> Iterable[float]: ...
    def __getitem__(self, key: int) -> float: ...


def vec(x: float, y: float, z: float) -> Vector: ...
def quat(x: float, y: float, z: float, w: float) -> Quaternion: ...
def mat(position: Vector = vec(0.0, 0.0, 0.0), rotation: Quaternion = quat(0.0, 0.0, 0.0, 1.0), scale: Vector = vec(1.0, 1.0, 1.0)) -> Matrix: ...

def rotate_x(angle: float) -> Quaternion: ...
def rotate_y(angle: float) -> Quaternion: ...
def rotate_z(angle: float) -> Quaternion: ...

def rotate(axis: Vector, angle: float) -> Quaternion: ...
def slerp(a: Quaternion, b: Quaternion, t: float) -> Quaternion: ...

def random_axis() -> Vector: ...
def random_rotation() -> Quaternion: ...