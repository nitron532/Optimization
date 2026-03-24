
Iteration has the form
$$x_{k+1} = x_k + \alpha_k p_k$$

And the search direction $p_k$ is required to be a descent direction for most line search algorithms, meaning that 
$$p_k^T\nabla f_k < 0$$
The search direction also often has the form
$$p_k^T = -B_k^{-1}\nabla f_k$$
where $B_k$ is a symmetric and non singular matrix. 
* In the steepest descent method, it is the identity matrix. 
* In Newton method, it is the exact Hessian. 
* In Quasi-Newton methods, it is an approximation of the Hessian that is updated at every iteration.


### Finding the ideal step length $\alpha$<hr>
Ideally, we would like to find the global minimizer of the univariate function 
$$\phi(\alpha) = f(x_k + \alpha p_k)$$
subject to $\alpha$ > 0.

It is often too expensive to find a local minimizer, let alone a global minimizer, which requires many expensive evaluations of the objective $f$ and its gradient.

Instead, an inexact line search is used to find a good enough $\alpha$. It is performed in two phases:
1. A bracketing phase, which finds an interval of desirable step lengths.
2. A bisection / interpolation phase, which finds a good step length in this interval

This is done until a condition is met, like an acceptable reduction in $f$.
### Only a reduction in $f$ is not enough to achieve convergence<hr>
Requiring only a reduction in $f$ is not enough to produce convergence to a minimum. Functions such as $f(x_k) = 5/k$ approach a limit at 0, but the minimum of this convex function is at -1:
![432](Images/Pasted%20image%2020260324133546.png)
To avoid this, we have to aim for a sufficient reduction in $f$, not just a reduction.

## Wolfe conditions <hr>

The Armijo condition says that $\alpha$ should satisfy this inequality:
$$f(x_k + \alpha p_k) \leq f(x_k) + c_1\alpha\nabla f_k^T p_k$$
where $c_1 \in (0,1)$ . 
Denote the right side of Armijo condition as $l(\alpha)$. It is the equation of the tangent of $f$ at point $x_k$along direction $p_k$. 
$c_1$ must not be 0 because then we only require the new iterate to be less than or equal to the previous, which clearly will not converge.
$c_1$ must also not be 1, because if the function is convex, the tangent is below the function, and the new iterate can not be under the tangent if the entire left hand side function is above the right hand side.
As $c_1 \rightarrow 0$, there are more available $\alpha$ to experiment with, and longer steps are allowed.
As $c_1 \rightarrow 1$, the interval is constricted to smaller steps, which are closer to the tangent (where there are no steps because the tangent is under a convex function, and the inequality can not be satisfied.)

It is a linear function with negative slope $c_1\nabla f_k^T p_k$ . The constant $c_1$, in practice chosen to be quite small (e.g $10^-4$) allows $l(\alpha) \geq \phi (\alpha)$ for small positive $\alpha$. This condition requires that the decrease in $f$ is proportional to the step length $\alpha$ and the directional derivative $\nabla f^T_k p_k$. Small $\alpha$ are expected to have small reductions, and large $\alpha$ are expected to have large reductions.

However, since the Armijo condition is valid for unreasonably small $\alpha$, a second condition known as the curvature condition, is imposed that constrains:
$$\nabla f(x_k+\alpha_k p_k)^Tp_k \geq c_2\nabla f_k^T p_k$$
The left hand side is simply $\phi'(\alpha)$ , saying that the slope at the new iterate should be greater than or equal to the previous iterate's slope. 
Geometrically, as we approach a minimum, the slope becomes more gentle, from a more steep negative, downhill, to a gentler, flatter, less negative slope (because we are reaching the bottom of the hill).
The parameter $c_2 \in (0,1)$.
As $c_2 \rightarrow 0$, we would like to be closer to the optimum.
As $c_2 \rightarrow 1$, we allow a softer increase in slope from the previous iterate.
For Newton/Quasi-Newton, $c_2 = 0.9$ is a common choice, whereas for nonlinear conjugate gradient, $c_2 = 0.1$ is a common choice.
![454](Images/Pasted%20image%2020260324133616.png)
Areas with steep downwards slope are not acceptable, since that means we could make a farther step length and achieve greater reductions.
Once we arrive at areas with gentler / uphill slopes, we want to stop, because we can no longer gain reductions going in this direction.

The Wolfe conditions' constants are constrained $0 < c_1 < c_2 < 1$.
The strong Wolfe conditions change the curvature condition to not allow the new derivative to be too positive.

Multiplying the objective function by a scalar or making an affine change to the variables does not change the Wolfe Conditions.


## Goldstein conditions <hr>
The Goldstein conditions also ensure that the step length achieves sufficient reductions and is not too short. It can be expressed in the following inequality
$$f(x_k) +(1-c)\alpha_k\nabla f_k^Tp_k \leq f(x_k + \alpha_kp_k) \leq f(x_k) +c\alpha_k\nabla f_k^Tp_k $$
Where $c \in (0, \frac{1}{2})$. Often used in Newton-type methods, but not Quasi-Newton.



### Backtracking line search <hr>
Since any $\alpha$ that is sufficiently small satisfies the Armijo condition, we introduced the curvature condition to encourage larger step sizes if possible.
Using backtracking, we start at a large step $\alpha$ which usually violates the Armijo condition (not enough reduction in $f$ for the step size $\alpha$), then reduce until it satisfies the condition. 
This guarantees that we don't suffer from the original problem with the Armijo condition, but also get a sufficiently large decrease in $f$.
This method of finding $\alpha$ only uses the Armijo condition (sufficient step size). Following the algorithm: ![626](Images/Pasted%20image%2020260324133641.png)
We choose a large initial $\alpha$, then iteratively decrease it by a possibly changing $\rho$ until the Armijo condition is satisfied. By starting large, we achieve the largest possible valid-according-to-Armijo-constraint step size (largest when reducing by $\rho$ each iteration), gaining a sufficient decrease in $f$ while also avoiding too small of a step size, which the curvature condition (not used in this algorithm) helps to prevent.

