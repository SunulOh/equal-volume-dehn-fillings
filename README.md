# equal-volume-dehn-fillings

This project explores Dehn fillings of a given hyperbolic 3-manifold that yield the same volume, but are not explained by identified symmetries in the manifold's volume formula.

---

### 📦 Prerequisites

This project is designed to run inside **SageMath** with the **SnapPy** package installed.

- SageMath version: `10.1`
- SnapPy version: `3.1.1`

🔗 Installation guides:
- [Install SageMath](https://doc.sagemath.org/html/en/installation/)
- [Install SnapPy](https://snappy.math.uic.edu/installing.html)

---

### ▶️ Usage Example

Once inside a SageMath environment, and with the repository cloned, try the following:

```python
import snappy
import symmetries
import volume

# Select first 10 manifolds from the cusped census
manifolds = snappy.OrientableCuspedCensus(num_cusps=1)[:10]

# Search for equal-volume Dehn fillings with coefficient range up to 100
volume.test(manifolds, 100)
```

To check whether specific Dehn filling coefficients preserve volume under the identified symmetries listed in `symmetries.txt`:

```python
symmetries.check(7, 11)
```
