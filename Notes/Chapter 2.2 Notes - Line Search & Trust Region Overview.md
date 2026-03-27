

Sylvester's Criterion:
* A way of proving the positive definitiveness of an n x n Hermitian matrix
* The upper left 1x1 determinant > 0, upper left 2x2 det. > 0, ... nxn det. >0

An iterate $x_k$ is the current estimation of the state.

Convex functions' local minima are also global minima. If there are multiple local minima in a convex function $x_1, x_2 ...$ all $f(x_1), f(x_2) ...$ are equal.

Note that $$\min_{p} f(x+p)$$ means that we minimize the objective $f$ by changing $p$.

## Line Search <hr>
Line search chooses a direction $p_k$ to step along, and the step length $\alpha$ from the current iterate $x_k$ can be found by approximating the one dimensional minimization $$\min_{\alpha > 0} f(x_k + \alpha p_k)$$
where $f$ is the objective function.
At $x_{k+1}$ a new search direction and step length is computed using this minimization. Often an exact solution is expensive and an approximation is enough. Usually a number of trial step lengths are attempted until one that loosely approximates the minimum is found.

Descent direction can be determined in many ways. If a direction $p_k$ is downhill, it satisfies $$p_k^T\nabla f_k = ||p_k||\ ||\nabla f_k||\cos(\theta_k) < 0$$ where $\theta_k$ is the angle between the descent direction $p_k$ and the gradient $\nabla f_k$.

### Gradient descent (steepest descent) <hr>
$p_k = -\nabla f$, the direction of steepest descent. Can be slow on more complicated problems

### Newton (second order descent direction) <hr>
Uses a quadratic model $$m_k(p) = f_k + p^T \nabla f_k + \frac{1}{2} p^T \nabla ^2f p \approx f(x+p)$$ then find $p$ that minimizes $m_k(p)$ using the explicit formula $$p_k^N = -(\nabla^2f_k)^{-1} \nabla f_k$$ 
	- If the Hessian is positive definite, the descent direction is then given by $\nabla f_k^T p^N_k$  as long as the gradient $\nabla f^T_k \neq 0$, which would also mean $p^N_k = 0$.  
	- The Hessian must be positive definite, or else $-(\nabla^2f_k)^{-1}$ will not exist, and the step direction is undefined. It also must be explicitly computed, which is a drawback for the Newton method.
Reliable when the difference between the objective $f(x_k+p)$ and its quadratic model $m_k(p)$ isn't too large.
Step size is usually $\alpha = 1$ and is only modified if insufficient results are produced

### Quasi-Newton <hr>
Quasi-Newton uses an approximation of the Hessian $B_k$ in exact place of the Hessian in the Newton method. So this gives:
$$p_k = -B_k^{-1} \nabla f_k$$
The next approximation $B_{k+1}$ must satisfy the secant equation $$B_{k+1}s_k = y_k $$where
- $s_k$ = $x_{k+1} - x_k$ 
- $y_k = \nabla f_{k+1} - \nabla f_k$

Other constraints may be imposed, such as symmetry, and that the difference between two approximations has low rank.
$B_{k+1}$ can be obtained from the following two formulas:
Symmetric rank-one formula, defined by
$$B_{k+1} = B_k + \frac{(y_k-B_ks_k)(y_k-B_ks_k)^T}{(y_k-B_ks_k)^Ts_k}$$
Or, the BFGS formula, defined by
$$B_{k+1} = B_k - \frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k} + \frac{y_ky_k^T}{y_k^Ts_k}$$
Both produce approximations that satisfy the secant equation.
### Conjugate gradient <hr>
[[Chapter 5 Notes]]



Any descent direction that makes an angle of $\leq \frac{\pi}{2}$ radians with $- \nabla f$ is guaranteed to decrease the objective function $f$ given that the step size is small enough.

## Trust Region <hr>
Trust region creates a model function $m_k$ that approximates the behavior of the actual objective function $f$ in the area near the current iterate $x_k$. 
$m_k$'s behavior may stray from the objective $f$'s the farther we go from $x_k$, so we restrict the search for the minimizer $x_k +p$ to a region around $x_k$, called the trust region.
The trust region radius $\Delta > 0$ constrains the step size to $||p||_2 \leq \Delta$  in the minimization $$\min_p m_k(x_k + p)$$If we find that the solution doesn't produce a good enough decrease in $f$, we reduce $\Delta$ and recompute the above equation.
This is because we assume that the behavior of $m_k$ near the edge of the trust region differed from the objective $f$ too greatly and so we must reduce $\Delta$ in order to restrain ourselves to a trust region that models $f$ accurately.

$m_k$ is usually defined as a quadratic of the form $$ m_k = f_k + p^T \nabla f_k + \frac{1}{2} p^T B_k p $$$f_k$ is the objective evaluated at the current iterate $x_k$
$\nabla f_k$ the gradient of the iterate at $x_k$
$B_k$ is either the Hessian at $x_k$ or some approximation of it.
These are chosen so that $m_k$ is in agreement with $f$ up to the second order **(?, textbook says only up to first order, is this because $B_k$ can be an approximation of the Hessian?)**

To determine the step direction multiple methods can be used:

### Line Search equivalent Method<hr>
If we set $B_k = 0$,  the minimization is now $$m_k(p) = f_k + p^T \nabla f_k$$
 subject to $||p||_2 \leq \Delta_k$. The closed form solution ends up being $$p_k = \frac{-\Delta_k \nabla f_k}{||\nabla f_k||} $$
which is equivalent to the line search gradient descent method, with a step size of the trust region radius $\Delta_k$.

## Scaling <hr>
A problem in unconstrained optimization is poorly scaled if changes to $x$ in one direction produce much larger variations in $f$ than they do with changes in other directions. For example, $$ f = 10^{100}x_2^2 + x_1^2 $$ is poorly scaled because a change in $x_2$ changes $f$ way more than a change in $x_1$ would.

Scaling can occur when we change the units of the variables, such as from meters to millimeters.
Steepest descent suffers when there is poor scaling, while Newton's method doesn't.

Generally speaking, incorporating scale invariance into line search algorithms is easier than for trust region algorithms.


