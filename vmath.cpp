#include <Python.h>
#include <structmember.h>

#include <random>
#include <chrono>

std::mt19937 mt((unsigned)std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> rng(0.0, 1.0);
const double pi = 3.1415926535897932;

struct Vector {
    PyObject_HEAD
    double x, y, z;
};

struct Quaternion {
    PyObject_HEAD
    double x, y, z, w;
};

struct Matrix {
    PyObject_HEAD
    double m[12];
};

PyTypeObject * Vector_type;
PyTypeObject * Quaternion_type;
PyTypeObject * Matrix_type;

bool is_float(PyObject * obj) {
    return Py_TYPE(obj)->tp_as_number && Py_TYPE(obj)->tp_as_number->nb_float;
}

double get_float(PyObject * obj) {
    if (PyFloat_CheckExact(obj)) {
        return PyFloat_AsDouble(obj);
    }
    if (PyLong_CheckExact(obj)) {
        return PyLong_AsDouble(obj);
    }
    PyObject * tmp = Py_TYPE(obj)->tp_as_number->nb_float(obj);
    double res = PyFloat_CheckExact(tmp) ? PyFloat_AsDouble(tmp) : nan(NULL);
    Py_DECREF(tmp);
    return res;
}

Vector * meth_vec(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"x", "y", "z", NULL};

    double x, y, z;

    if (PyTuple_Size(args) == 1 && (!kwargs || PyDict_Size(kwargs) == 0)) {
        if (!PyArg_ParseTuple(args, "(ddd)", &x, &y, &z)) {
            return NULL;
        }
    } else if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddd", keywords, &x, &y, &z)) {
        return NULL;
    }

    Vector * res = PyObject_New(Vector, Vector_type);
    res->x = x;
    res->y = y;
    res->z = z;
    return res;
}

Quaternion * meth_quat(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"x", "y", "z", "w", NULL};

    double x, y, z, w;

    if (PyTuple_Size(args) == 1 && (!kwargs || PyDict_Size(kwargs) == 0)) {
        if (!PyArg_ParseTuple(args, "(dddd)", &x, &y, &z, &w)) {
            return NULL;
        }
    } else if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddd", keywords, &x, &y, &z, &w)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double lng = sqrt(x * x + y * y + z * z + w * w) * (w < 0.0 ? -1.0 : 1.0);
    res->x = x / lng;
    res->y = y / lng;
    res->z = z / lng;
    res->w = w / lng;
    return res;
}

Matrix * meth_mat(PyObject * self, PyObject * args, PyObject * kwargs) {
    if (PyTuple_Size(args) == 1 && (!kwargs || PyDict_Size(kwargs) == 0)) {
        double m[12];
        if (!PyArg_ParseTuple(args, "(dddddddddddd)", &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8], &m[9], &m[10], &m[11])) {
            return NULL;
        }
        Matrix * res = PyObject_New(Matrix, Matrix_type);
        memcpy(res->m, m, sizeof(m));
        return res;
    }

    static char * keywords[] = {"position", "rotation", "scale", NULL};

    static Vector default_position = {{}, 0.0, 0.0, 0.0};
    static Quaternion default_rotation = {{}, 0.0, 0.0, 0.0, 1.0};
    static Vector default_scale = {{}, 1.0, 1.0, 1.0};

    Vector * p = &default_position;
    Quaternion * r = &default_rotation;
    Vector * s = &default_scale;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|O!O!O!", keywords, Vector_type, &p, Quaternion_type, &r, Vector_type, &s)) {
        return NULL;
    }

    Matrix * res = PyObject_New(Matrix, Matrix_type);
    res->m[0] = (1.0 - 2.0 * r->y * r->y - 2.0 * r->z * r->z) * s->x;
    res->m[1] = (2.0 * r->x * r->y - 2.0 * r->z * r->w) * s->x;
    res->m[2] = (2.0 * r->x * r->z + 2.0 * r->y * r->w) * s->x;
    res->m[3] = p->x;
    res->m[4] = (2.0 * r->x * r->y + 2.0 * r->z * r->w) * s->y;
    res->m[5] = (1.0 - 2.0 * r->x * r->x - 2.0 * r->z * r->z) * s->y;
    res->m[6] = (2.0 * r->y * r->z - 2.0 * r->x * r->w) * s->y;
    res->m[7] = p->y;
    res->m[8] = (2.0 * r->x * r->z - 2.0 * r->y * r->w) * s->z;
    res->m[9] = (2.0 * r->y * r->z + 2.0 * r->x * r->w) * s->z;
    res->m[10] = (1.0 - 2.0 * r->x * r->x - 2.0 * r->y * r->y) * s->z;
    res->m[11] = p->z;
    return res;
}

Vector * meth_random_axis(PyObject * self) {
    Vector * res = PyObject_New(Vector, Vector_type);
    while (true) {
        const double u1 = rng(mt) * 2.0 - 1.0;
        const double u2 = rng(mt) * 2.0 - 1.0;
        const double u3 = rng(mt) * 2.0 - 1.0;
        if (u1 * u1 + u2 * u2 + u3 * u3 > 1.0) {
            continue;
        }
        res->x = u1;
        res->y = u2;
        res->z = u3;
        break;
    }
    return res;
}

Quaternion * meth_random_rotation(PyObject * self) {
    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double u1 = rng(mt);
    const double u2 = rng(mt);
    const double u3 = rng(mt);
    res->x = sqrt(1.0 - u1) * sin(2.0 * pi * u2);
    res->y = sqrt(1.0 - u1) * cos(2.0 * pi * u2);
    res->z = sqrt(u1) * sin(2.0 * pi * u3);
    res->w = sqrt(u1) * cos(2.0 * pi * u3);
    if (res->w < 0.0) {
        res->x = -res->x;
        res->y = -res->y;
        res->z = -res->z;
        res->w = -res->w;
    }
    return res;
}

Quaternion * meth_rotate_x(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"angle", NULL};

    double angle;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", keywords, &angle)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double s = sin(angle / 2.0);
    const double c = cos(angle / 2.0);
    const double sgn = c < 0.0 ? -1.0 : 1.0;
    res->x = s * sgn;
    res->y = 0.0;
    res->z = 0.0;
    res->w = c * sgn;
    return res;
}

Quaternion * meth_rotate_y(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"angle", NULL};

    double angle;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", keywords, &angle)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double s = sin(angle / 2.0);
    const double c = cos(angle / 2.0);
    const double sgn = c < 0.0 ? -1.0 : 1.0;
    res->x = 0.0;
    res->y = s * sgn;
    res->z = 0.0;
    res->w = c * sgn;
    return res;
}

Quaternion * meth_rotate_z(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"angle", NULL};

    double angle;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", keywords, &angle)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double s = sin(angle / 2.0);
    const double c = cos(angle / 2.0);
    const double sgn = c < 0.0 ? -1.0 : 1.0;
    res->x = 0.0;
    res->y = 0.0;
    res->z = s * sgn;
    res->w = c * sgn;
    return res;
}

Quaternion * Quaternion_meth_inverse(Quaternion * self) {
    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    res->x = -self->x;
    res->y = -self->y;
    res->z = -self->z;
    res->w = self->w;
    return res;
}

Matrix * Matrix_meth_inverse(Matrix * self) {
    const double t1 = self->m[5] * self->m[10] - self->m[9] * self->m[6];
    const double t2 = self->m[9] * self->m[2] - self->m[1] * self->m[10];
    const double t3 = self->m[1] * self->m[6] - self->m[5] * self->m[2];
    const double d = 1.0 / (self->m[0] * t1 + self->m[4] * t2 + self->m[8] * t3);

    Matrix * res = PyObject_New(Matrix, Matrix_type);
    res->m[0] = t1 * d;
    res->m[4] = (self->m[8] * self->m[6] - self->m[4] * self->m[10]) * d;
    res->m[8] = (self->m[4] * self->m[9] - self->m[8] * self->m[5]) * d;
    res->m[1] = t2 * d;
    res->m[5] = (self->m[0] * self->m[10] - self->m[8] * self->m[2]) * d;
    res->m[9] = (self->m[8] * self->m[1] - self->m[0] * self->m[9]) * d;
    res->m[2] = t3 * d;
    res->m[6] = (self->m[4] * self->m[2] - self->m[0] * self->m[6]) * d;
    res->m[10] = (self->m[0] * self->m[5] - self->m[4] * self->m[1]) * d;
    res->m[3] = -(res->m[0] * self->m[3] + res->m[1] * self->m[7] + res->m[2] * self->m[11]);
    res->m[7] = -(res->m[4] * self->m[3] + res->m[5] * self->m[7] + res->m[6] * self->m[11]);
    res->m[11] = -(res->m[8] * self->m[3] + res->m[9] * self->m[7] + res->m[10] * self->m[11]);
    return res;
}

Vector * Matrix_meth_position(Matrix * self) {
    Vector * res = PyObject_New(Vector, Vector_type);
    res->x = self->m[3];
    res->y = self->m[7];
    res->z = self->m[11];
    return res;
}

Quaternion * Matrix_meth_rotation(Matrix * self) {
    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double cx = self->m[1] * self->m[6] - self->m[5] * self->m[2];
    const double cy = self->m[2] * self->m[4] - self->m[6] * self->m[0];
    const double cz = self->m[0] * self->m[5] - self->m[4] * self->m[1];
    const double bx = cy * self->m[2] - self->m[1] * cz;
    const double by = cz * self->m[0] - self->m[2] * cx;
    const double bz = cx * self->m[1] - self->m[0] * cy;
    const double la = sqrt(self->m[0] * self->m[0] + self->m[1] * self->m[1] + self->m[2] * self->m[2]);
    const double lb = sqrt(bx * bx + by * by + bz * bz);
    const double lc = sqrt(cx * cx + cy * cy + cz * cz);
    const double m[9] = {
        self->m[0] / la, self->m[1] / la, self->m[2] / la,
        bx / lb, by / lb, bz / lb,
        cx / lc, cy / lc, cz / lc,
    };
    const double sx = m[7] < m[5] ? -1.0 : 1.0;
    const double sy = m[2] < m[6] ? -1.0 : 1.0;
    const double sz = m[3] < m[1] ? -1.0 : 1.0;
    const double tx = 1.0 + m[0] - m[4] - m[8];
    const double ty = 1.0 - m[0] + m[4] - m[8];
    const double tz = 1.0 - m[0] - m[4] + m[8];
    const double tw = 1.0 + m[0] + m[4] + m[8];
    res->x = tx > 0.0 ? sqrt(tx) / 2.0 * sx : 0.0;
    res->y = ty > 0.0 ? sqrt(ty) / 2.0 * sy : 0.0;
    res->z = tz > 0.0 ? sqrt(tz) / 2.0 * sz : 0.0;
    res->w = tw > 0.0 ? sqrt(tw) / 2.0 : 0.0;
    return res;
}

Vector * Matrix_meth_scale(Matrix * self) {
    Vector * res = PyObject_New(Vector, Vector_type);
    res->x = sqrt(self->m[0] * self->m[0] + self->m[1] * self->m[1] + self->m[2] * self->m[2]);
    res->y = sqrt(self->m[4] * self->m[4] + self->m[5] * self->m[5] + self->m[6] * self->m[6]);
    res->z = sqrt(self->m[8] * self->m[8] + self->m[9] * self->m[9] + self->m[10] * self->m[10]);
    return res;
}

Vector * Quaternion_meth_axis(Quaternion * self) {
    Vector * res = PyObject_New(Vector, Vector_type);
    const double t = sqrt(1.0 - self->w * self->w);
    if (t) {
        res->x = self->x / t;
        res->y = self->y / t;
        res->z = self->z / t;
    } else {
        res->x = 0.0;
        res->y = 0.0;
        res->z = 1.0;
    }
    return res;
}

PyObject * Quaternion_meth_angle(Quaternion * arg) {
    return PyFloat_FromDouble(acos(arg->w) * 2.0);
}

Vector * Vector_meth_normal(Vector * self) {
    const double lng = sqrt(self->x * self->x + self->y * self->y + self->z * self->z);
    Vector * res = PyObject_New(Vector, Vector_type);
    if (lng) {
        res->x = self->x / lng;
        res->y = self->y / lng;
        res->z = self->z / lng;
    } else {
        res->x = 0.0;
        res->y = 0.0;
        res->z = 1.0;
    }
    return res;
}

PyObject * Vector_meth_length(Vector * self) {
    return PyFloat_FromDouble(sqrt(self->x * self->x + self->y * self->y + self->z * self->z));
}

Quaternion * meth_rotate(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"axis", "angle", NULL};

    Vector * axis;
    double angle;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!d", keywords, Vector_type, &axis, &angle)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
    const double s = sin(angle / 2.0) / sqrt(axis->x * axis->x + axis->y * axis->y + axis->z * axis->z);
    const double x = axis->x * s;
    const double y = axis->y * s;
    const double z = axis->z * s;
    const double w = cos(angle / 2.0);
    const double sgn = w < 0.0 ? -1.0 : 1.0;
    res->x = x * sgn;
    res->y = y * sgn;
    res->z = z * sgn;
    res->w = w * sgn;
    return res;
}

Quaternion * meth_slerp(PyObject * self, PyObject * args, PyObject * kwargs) {
    static char * keywords[] = {"a", "b", "t", NULL};

    Quaternion * a;
    Quaternion * b;
    double t;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!d", keywords, Quaternion_type, &a, Quaternion_type, &b, &t)) {
        return NULL;
    }

    Quaternion * res = PyObject_New(Quaternion, Quaternion_type);

    const double qx = a->w * b->x - a->x * b->w - a->y * b->z + a->z * b->y;
    const double qy = a->w * b->y - a->y * b->w - a->z * b->x + a->x * b->z;
    const double qz = a->w * b->z - a->z * b->w - a->x * b->y + a->y * b->x;
    const double qw = a->w * b->w + a->x * b->x + a->y * b->y + a->z * b->z;
    const double angle = acos(qw) * 2.0 * t;

    if (!angle) {
        res->x = a->x;
        res->y = a->y;
        res->z = a->z;
        res->w = a->w;
        return res;
    }

    const double p = sqrt(1.0 - qw * qw);
    const double tx = qx / p;
    const double ty = qy / p;
    const double tz = qz / p;

    const double s = sin(angle / 2.0) / sqrt(tx * tx + ty * ty + tz * tz);
    const double rx = tx * s;
    const double ry = ty * s;
    const double rz = tz * s;
    const double rw = cos(angle / 2.0);

    res->x = rw * a->x + rx * a->w + ry * a->z - rz * a->y;
    res->y = rw * a->y + ry * a->w + rz * a->x - rx * a->z;
    res->z = rw * a->z + rz * a->w + rx * a->y - ry * a->x;
    res->w = rw * a->w - rx * a->x - ry * a->y - rz * a->z;

    const double lng = sqrt(res->x * res->x + res->y * res->y + res->z * res->z + res->w * res->w) * (res->w < 0.0 ? -1.0 : 1.0);
    res->x = res->x / lng;
    res->y = res->y / lng;
    res->z = res->z / lng;
    res->w = res->w / lng;
    return res;
}

PyObject * Vector_nb_add(Vector * self, PyObject * other) {
    if (Py_TYPE(self) == Vector_type && Py_TYPE(other) == Vector_type) {
        Vector * b = (Vector *)other;
        Vector * res = PyObject_New(Vector, Vector_type);
        res->x = self->x + b->x;
        res->y = self->y + b->y;
        res->z = self->z + b->z;
        return (PyObject *)res;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Vector_nb_subtract(Vector * self, PyObject * other) {
    if (Py_TYPE(self) == Vector_type && Py_TYPE(other) == Vector_type) {
        Vector * b = (Vector *)other;
        Vector * res = PyObject_New(Vector, Vector_type);
        res->x = self->x - b->x;
        res->y = self->y - b->y;
        res->z = self->z - b->z;
        return (PyObject *)res;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Vector_nb_negative(Vector * self) {
    Vector * res = PyObject_New(Vector, Vector_type);
    res->x = -self->x;
    res->y = -self->y;
    res->z = -self->z;
    return (PyObject *)res;
}

PyObject * Vector_nb_multiply(PyObject * self, PyObject * other) {
    if (Py_TYPE(self) == Vector_type) {
        if (Py_TYPE(other) == Vector_type) {
            Vector * a = (Vector *)self;
            Vector * b = (Vector *)other;
            Vector * res = PyObject_New(Vector, Vector_type);
            res->x = a->x * b->x;
            res->y = a->y * b->y;
            res->z = a->z * b->z;
            return (PyObject *)res;
        }
        if (is_float(other)) {
            Vector * a = (Vector *)self;
            const double b = get_float(other);
            Vector * res = PyObject_New(Vector, Vector_type);
            res->x = a->x * b;
            res->y = a->y * b;
            res->z = a->z * b;
            return (PyObject *)res;
        }
    }
    if (is_float(self) && Py_TYPE(other) == Vector_type) {
        const double a = get_float(other);
        Vector * b = (Vector *)self;
        Vector * res = PyObject_New(Vector, Vector_type);
        res->x = a * b->x;
        res->y = a * b->y;
        res->z = a * b->z;
        return (PyObject *)res;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Vector_nb_true_divide(PyObject * self, PyObject * other) {
    if (Py_TYPE(self) == Vector_type) {
        if (Py_TYPE(other) == Vector_type) {
            Vector * a = (Vector *)self;
            Vector * b = (Vector *)other;
            Vector * res = PyObject_New(Vector, Vector_type);
            res->x = a->x / b->x;
            res->y = a->y / b->y;
            res->z = a->z / b->z;
            return (PyObject *)res;
        }
        if (is_float(other)) {
            Vector * a = (Vector *)self;
            const double b = get_float(other);
            Vector * res = PyObject_New(Vector, Vector_type);
            res->x = a->x / b;
            res->y = a->y / b;
            res->z = a->z / b;
            return (PyObject *)res;
        }
    }
    if (is_float(self) && Py_TYPE(other) == Vector_type) {
        const double a = get_float(self);
        Vector * b = (Vector *)other;
        Vector * res = PyObject_New(Vector, Vector_type);
        res->x = a / b->x;
        res->y = a / b->y;
        res->z = a / b->z;
        return (PyObject *)res;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Quaternion_nb_multiply(PyObject * self, PyObject * other) {
    if (Py_TYPE(self) == Quaternion_type) {
        if (Py_TYPE(other) == Quaternion_type) {
            Quaternion * a = (Quaternion *)self;
            Quaternion * b = (Quaternion *)other;
            Quaternion * res = PyObject_New(Quaternion, Quaternion_type);
            res->x = a->w * b->x + a->x * b->w + a->y * b->z - a->z * b->y;
            res->y = a->w * b->y + a->y * b->w + a->z * b->x - a->x * b->z;
            res->z = a->w * b->z + a->z * b->w + a->x * b->y - a->y * b->x;
            res->w = a->w * b->w - a->x * b->x - a->y * b->y - a->z * b->z;
            const double lng = sqrt(res->x * res->x + res->y * res->y + res->z * res->z + res->w * res->w) * (res->w < 0.0 ? -1.0 : 1.0);
            res->x = res->x / lng;
            res->y = res->y / lng;
            res->z = res->z / lng;
            res->w = res->w / lng;
            return (PyObject *)res;
        }
        if (Py_TYPE(other) == Vector_type) {
            Quaternion * a = (Quaternion *)self;
            Vector * b = (Vector *)other;
            Vector * res = PyObject_New(Vector, Vector_type);
            const double tx = b->y * a->z - a->y * b->z - a->w * b->x;
            const double ty = b->z * a->x - a->z * b->x - a->w * b->y;
            const double tz = b->x * a->y - a->x * b->y - a->w * b->z;
            res->x = b->x + (ty * a->z - a->y * tz) * 2.0;
            res->y = b->y + (tz * a->x - a->z * tx) * 2.0;
            res->z = b->z + (tx * a->y - a->x * ty) * 2.0;
            return (PyObject *)res;
        }
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Matrix_nb_multiply(PyObject * self, PyObject * other) {
    if (Py_TYPE(self) == Matrix_type) {
        if (Py_TYPE(other) == Matrix_type) {
            Matrix * a = (Matrix *)self;
            Matrix * b = (Matrix *)other;
            Matrix * res = PyObject_New(Matrix, Matrix_type);
            res->m[0] = a->m[0] * b->m[0] + a->m[1] * b->m[4] + a->m[2] * b->m[8];
            res->m[1] = a->m[0] * b->m[1] + a->m[1] * b->m[5] + a->m[2] * b->m[9];
            res->m[2] = a->m[0] * b->m[2] + a->m[1] * b->m[6] + a->m[2] * b->m[10];
            res->m[3] = a->m[0] * b->m[3] + a->m[1] * b->m[7] + a->m[2] * b->m[11] + a->m[3];
            res->m[4] = a->m[4] * b->m[0] + a->m[5] * b->m[4] + a->m[6] * b->m[8];
            res->m[5] = a->m[4] * b->m[1] + a->m[5] * b->m[5] + a->m[6] * b->m[9];
            res->m[6] = a->m[4] * b->m[2] + a->m[5] * b->m[6] + a->m[6] * b->m[10];
            res->m[7] = a->m[4] * b->m[3] + a->m[5] * b->m[7] + a->m[6] * b->m[11] + a->m[7];
            res->m[8] = a->m[8] * b->m[0] + a->m[9] * b->m[4] + a->m[10] * b->m[8];
            res->m[9] = a->m[8] * b->m[1] + a->m[9] * b->m[5] + a->m[10] * b->m[9];
            res->m[10] = a->m[8] * b->m[2] + a->m[9] * b->m[6] + a->m[10] * b->m[10];
            res->m[11] = a->m[8] * b->m[3] + a->m[9] * b->m[7] + a->m[10] * b->m[11] + a->m[11];
            return (PyObject *)res;
        }
        if (Py_TYPE(other) == Vector_type) {
            Matrix * a = (Matrix *)self;
            Vector * b = (Vector *)other;
            Vector * res = PyObject_New(Vector, Vector_type);
            res->x = a->m[0] * b->x + a->m[1] * b->y + a->m[2] * b->z + a->m[3];
            res->y = a->m[4] * b->x + a->m[5] * b->y + a->m[6] * b->z + a->m[7];
            res->z = a->m[8] * b->x + a->m[9] * b->y + a->m[10] * b->z + a->m[11];
            return (PyObject *)res;
        }
    }
    Py_RETURN_NOTIMPLEMENTED;
}

PyObject * Vector_meth_pack(Vector * self) {
    PyObject * res = PyBytes_FromStringAndSize(NULL, 12);
    float * ptr = (float *)PyBytes_AsString(res);
    ptr[0] = (float)self->x;
    ptr[1] = (float)self->y;
    ptr[2] = (float)self->z;
    return res;
}

PyObject * Quaternion_meth_pack(Quaternion * self) {
    PyObject * res = PyBytes_FromStringAndSize(NULL, 16);
    float * ptr = (float *)PyBytes_AsString(res);
    ptr[0] = (float)self->x;
    ptr[1] = (float)self->y;
    ptr[2] = (float)self->z;
    ptr[3] = (float)self->w;
    return res;
}

PyObject * Matrix_meth_pack(Matrix * self) {
    PyObject * res = PyBytes_FromStringAndSize(NULL, 64);
    float * ptr = (float *)PyBytes_AsString(res);
    for (int i = 0; i < 12; ++i) {
        ptr[i] = (float)self->m[i];
    }
    ptr[12] = 0.0f;
    ptr[13] = 0.0f;
    ptr[14] = 0.0f;
    ptr[15] = 1.0f;
    return res;
}

Py_ssize_t Vector_sq_length(Vector * self) {
    return 3;
}

Py_ssize_t Quaternion_sq_length(Quaternion * self) {
    return 4;
}

Py_ssize_t Matrix_sq_length(Matrix * self) {
    return 12;
}

PyObject * Vector_sq_item(Vector * self, Py_ssize_t idx) {
    switch (idx) {
        case 0: return PyFloat_FromDouble(self->x);
        case 1: return PyFloat_FromDouble(self->y);
        case 2: return PyFloat_FromDouble(self->z);
    }
    PyErr_Format(PyExc_IndexError, "");
    return NULL;
}

PyObject * Quaternion_sq_item(Quaternion * self, Py_ssize_t idx) {
    switch (idx) {
        case 0: return PyFloat_FromDouble(self->x);
        case 1: return PyFloat_FromDouble(self->y);
        case 2: return PyFloat_FromDouble(self->z);
        case 3: return PyFloat_FromDouble(self->w);
    }
    PyErr_Format(PyExc_IndexError, "");
    return NULL;
}

PyObject * Matrix_sq_item(Matrix * self, Py_ssize_t idx) {
    if (idx >= 0 && idx < 12) {
        return PyFloat_FromDouble(self->m[idx]);
    }
    PyErr_Format(PyExc_IndexError, "");
    return NULL;
}

PyObject * default_str(PyObject * self) {
    PyObject * tmp = PySequence_Tuple(self);
    PyObject * res = PyObject_Str(tmp);
    Py_DECREF(tmp);
    return res;
}

void default_dealloc(PyObject * self) {
    Py_TYPE(self)->tp_free(self);
}

PyMethodDef Vector_methods[] = {
    {"normal", (PyCFunction)Vector_meth_normal, METH_NOARGS, NULL},
    {"length", (PyCFunction)Vector_meth_length, METH_NOARGS, NULL},
    {"pack", (PyCFunction)Vector_meth_pack, METH_NOARGS, NULL},
    {},
};

PyMethodDef Quaternion_methods[] = {
    {"axis", (PyCFunction)Quaternion_meth_axis, METH_NOARGS, NULL},
    {"angle", (PyCFunction)Quaternion_meth_angle, METH_NOARGS, NULL},
    {"inverse", (PyCFunction)Quaternion_meth_inverse, METH_NOARGS, NULL},
    {"pack", (PyCFunction)Quaternion_meth_pack, METH_NOARGS, NULL},
    {},
};

PyMethodDef Matrix_methods[] = {
    {"inverse", (PyCFunction)Matrix_meth_inverse, METH_NOARGS, NULL},
    {"position", (PyCFunction)Matrix_meth_position, METH_NOARGS, NULL},
    {"rotation", (PyCFunction)Matrix_meth_rotation, METH_NOARGS, NULL},
    {"scale", (PyCFunction)Matrix_meth_scale, METH_NOARGS, NULL},
    {"pack", (PyCFunction)Matrix_meth_pack, METH_NOARGS, NULL},
    {},
};

PyMemberDef Vector_members[] = {
    {"x", T_DOUBLE, offsetof(Vector, x), READONLY, NULL},
    {"y", T_DOUBLE, offsetof(Vector, y), READONLY, NULL},
    {"z", T_DOUBLE, offsetof(Vector, z), READONLY, NULL},
    {},
};

PyMemberDef Quaternion_members[] = {
    {"x", T_DOUBLE, offsetof(Quaternion, x), READONLY, NULL},
    {"y", T_DOUBLE, offsetof(Quaternion, y), READONLY, NULL},
    {"z", T_DOUBLE, offsetof(Quaternion, z), READONLY, NULL},
    {"w", T_DOUBLE, offsetof(Quaternion, w), READONLY, NULL},
    {},
};

PyType_Slot Vector_slots[] = {
    {Py_tp_methods, Vector_methods},
    {Py_tp_members, Vector_members},
    {Py_nb_add, Vector_nb_add},
    {Py_nb_subtract, Vector_nb_subtract},
    {Py_nb_negative, Vector_nb_negative},
    {Py_nb_multiply, Vector_nb_multiply},
    {Py_nb_true_divide, Vector_nb_true_divide},
    {Py_sq_length, Vector_sq_length},
    {Py_sq_item, Vector_sq_item},
    {Py_tp_dealloc, default_dealloc},
    {Py_tp_str, default_str},
    {},
};

PyType_Slot Quaternion_slots[] = {
    {Py_tp_methods, Quaternion_methods},
    {Py_tp_members, Quaternion_members},
    {Py_nb_multiply, Quaternion_nb_multiply},
    {Py_sq_length, Quaternion_sq_length},
    {Py_sq_item, Quaternion_sq_item},
    {Py_tp_dealloc, default_dealloc},
    {Py_tp_str, default_str},
    {},
};

PyType_Slot Matrix_slots[] = {
    {Py_tp_methods, Matrix_methods},
    {Py_nb_multiply, Matrix_nb_multiply},
    {Py_sq_length, Matrix_sq_length},
    {Py_sq_item, Matrix_sq_item},
    {Py_tp_dealloc, default_dealloc},
    {Py_tp_str, default_str},
    {},
};

PyType_Spec Vector_spec = {"vmath.Vector", sizeof(Vector), 0, Py_TPFLAGS_DEFAULT, Vector_slots};
PyType_Spec Quaternion_spec = {"vmath.Quaternion", sizeof(Quaternion), 0, Py_TPFLAGS_DEFAULT, Quaternion_slots};
PyType_Spec Matrix_spec = {"vmath.Matrix", sizeof(Matrix), 0, Py_TPFLAGS_DEFAULT, Matrix_slots};

PyMethodDef module_methods[] = {
    {"vec", (PyCFunction)meth_vec, METH_VARARGS | METH_KEYWORDS, NULL},
    {"quat", (PyCFunction)meth_quat, METH_VARARGS | METH_KEYWORDS, NULL},
    {"mat", (PyCFunction)meth_mat, METH_VARARGS | METH_KEYWORDS, NULL},
    {"rotate", (PyCFunction)meth_rotate, METH_VARARGS | METH_KEYWORDS, NULL},
    {"slerp", (PyCFunction)meth_slerp, METH_VARARGS | METH_KEYWORDS, NULL},
    {"random_axis", (PyCFunction)meth_random_axis, METH_NOARGS, NULL},
    {"random_rotation", (PyCFunction)meth_random_rotation, METH_NOARGS, NULL},
    {"rotate_x", (PyCFunction)meth_rotate_x, METH_VARARGS | METH_KEYWORDS, NULL},
    {"rotate_y", (PyCFunction)meth_rotate_y, METH_VARARGS | METH_KEYWORDS, NULL},
    {"rotate_z", (PyCFunction)meth_rotate_z, METH_VARARGS | METH_KEYWORDS, NULL},
    {},
};

PyModuleDef module_def = {PyModuleDef_HEAD_INIT, "vmath", NULL, -1, module_methods};

extern "C" PyObject * PyInit_vmath() {
    PyObject * module = PyModule_Create(&module_def);
    Vector_type = (PyTypeObject *)PyType_FromSpec(&Vector_spec);
    Quaternion_type = (PyTypeObject *)PyType_FromSpec(&Quaternion_spec);
    Matrix_type = (PyTypeObject *)PyType_FromSpec(&Matrix_spec);
    PyModule_AddObject(module, "pi", PyFloat_FromDouble(pi));
    return module;
}
