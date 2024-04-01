# Bregman first-order methods
This repository hosts an implementation of first-order methods for optimization over the cone of nonnegative trigonometric matrix polynomials $K$ given by

$$
K = \\{(X_0, X_1, \dots, X_p) \\: | \\:  X_0 + \sum_{k=1}^p (X_k z^k + X_k z^{-k}) \succeq 0, \\: \forall z \in \mathbb{C} : |z| = 1\\},
$$

where the matrices $X_0 \in \mathbb{S}^m$ and $X_k \in \mathbb{R}^{m \times m}$, $1 \leq k \leq p$ are variables in the optimization problem. For details and applications,
we refer to our manuscript XXX.

## Guide to code
There are essentially three main folders:
1. The Bregman proximal operator and the two algorithms in the paper referred to as Bregman Douglas--Rachford and Improved interior gradient algorithm are implemented in the `algorithms` folder. 
2. The results presented in Section 5.1 of the paper were produced using the code in the  `eucl_proj` folder.
3. The results presented in Section 5.2 of the paper were produced using the code in the  `graphical_models` folder.

If you have any questions, feel free to reach out.

## Citing
If you wish to cite this work, please use the following BibTeX:

```
@ARTICLE{CED24,
  author={Cederberg, D.},
  journal={TODO}, 
  title={First-order methods for optimization over nonnegative trigonometric matrix polynomials}, 
  year={2024}
}
```
