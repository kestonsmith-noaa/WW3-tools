import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, rand

# Example: Create a large sparse matrix Q and vector r
n = 100
Q = rand(n, n, density=0.01, format='csr') # CSR format is efficient for products
#Q = Q + Q.T # Make it symmetric
Q = Q @ Q.T # Make it symmetric
r = np.random.rand(n)

# Initial guess
x0 = np.zeros(n)
def objective_function(x):
    # Calculate 0.5 * x.T @ Q @ x + r.T @ x
    return 0.5 * x.T @ Q @ x + r.T @ x

def gradient_function(x):
    # Calculate Q @ x + r
    return Q @ x + r

def hessian_product_function(p, x_ignored):
    # Calculate H @ p, where H is Q. The 'x_ignored' argument is required by the API.
    return Q @ p

x0=x0+0.1
res = minimize(
    fun=objective_function, 
    x0=x0, 
    method='Newton-CG', 
    jac=gradient_function, 
    hessp=lambda p, x: Q @ p, # Or use the defined hessian_product_function
    options={'disp': True}
)

print("Optimization successful:", res.success)
print("Optimal x:", res.x)
print("Minimum value:", res.fun)
