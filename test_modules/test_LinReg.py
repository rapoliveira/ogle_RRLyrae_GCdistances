import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression 
from sklearn.metrics import mean_squared_error, r2_score

### Example from Geron's ML 2019 book (page 170)

X = 2 * np.random.rand(100,1)
y = 4 + 3*X + np.random.randn(100,1)
X_b = np.c_[np.ones((100,1)),X]     # add x0=1 to each instance

#.with NumPy's linear algebra 
theta_best = np.linalg.inv(X_b.T.dot(X_b)).dot(X_b.T).dot(y)

#.with Scikit-Learn 
lin_reg = LinearRegression()
lin_reg.fit(X, y)
print ('- Intercep =',lin_reg.intercept_[0])
print ('- Angular coef =', lin_reg.coef_[0][0])

X_new2 = np.array([[0],[2]])
aux = np.arange(0,2,0.02)
X_new = np.array([[aux[i]] for i in range(len(aux))])

X_new_b = np.c_[np.ones((100,1)),X_new]
y_predict = X_new_b.dot(theta_best)

print (y.shape,lin_reg.predict(X_new).shape)
input()
print('Mean squared error: %.2f'
      % mean_squared_error(y, lin_reg.predict(X_new)))
print('Coefficient of determination: %.2f'
      % r2_score(y, lin_reg.predict(X_new)))

fig = plt.figure(figsize=(6,4))
plt.plot(X_new, y_predict, 'r--', lw=5, label='NumPy')
plt.plot(X_new, lin_reg.predict(X_new), 'b--', lw=2, label='Scikit')
plt.scatter(X, y, color='gray',s=10)
plt.axis([0, 2, 0, 15])
plt.xlabel('x',fontsize=11)
plt.ylabel('y',fontsize=11)
plt.legend(loc=2,fontsize=11)
plt.tight_layout()
plt.show()


### Example from sklearn tutorial

'''# Load the diabetes dataset
diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

# Use only one feature
diabetes_X = diabetes_X[:, np.newaxis, 2]

print (diabetes_X)
print (len(diabetes_X),type(diabetes_X),diabetes_X.shape)

# Split the data into training/testing sets
diabetes_X_train = diabetes_X[:-20]
diabetes_X_test = diabetes_X[-20:]
print (len(diabetes_X_train),len(diabetes_X_test))

# Split the targets into training/testing sets
diabetes_y_train = diabetes_y[:-20]
diabetes_y_test = diabetes_y[-20:]

regr = linear_model.LinearRegression()
regr.fit(diabetes_X_train, diabetes_y_train)
diabetes_y_pred = regr.predict(diabetes_X_test)

print('Coefficients: \n', regr.coef_, regr.intercept_)
print('Mean squared error: %.2f'
      % mean_squared_error(diabetes_y_test, diabetes_y_pred))
print('Coefficient of determination: %.2f'
      % r2_score(diabetes_y_test, diabetes_y_pred))

plt.scatter(diabetes_X_test, diabetes_y_test,  color='black')
plt.plot(diabetes_X_test, diabetes_y_pred, color='blue', linewidth=3)

print (len(diabetes_X_test))

plt.tight_layout()
plt.show()'''