import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

x = np.linspace(-5, 5, 20)

x1 = x*np.random.rand(20)
x2 = x*np.random.rand(20) + 10
x3 = x*np.random.rand(20) - 10

xt = np.hstack((x1,x2,x3))  # horizontal combination
plt.hist(xt)
plt.show()

max_iterations = 10
pi = np.array([1/3, 1/3, 1/3])
mu = np.array([5,6,-3])
var = np.array([1,3,9])
r = np.zeros((len(xt), 3))

gauss1 = norm(loc=mu[0], scale=var[0])
gauss2 = norm(loc=mu[1], scale=var[1])
gauss3 = norm(loc=mu[2], scale=var[2])

# E-Step
for c,g,p in zip(range(3), [gauss1, gauss2, gauss3], pi):
    r[:,c] = p*g.pdf(xt[:])

for i in range(len(r)):
    r[i,:] /= np.sum(r[i,:])

fig = plt.figure(figsize=(6,5))
ax0 = fig.add_subplot(111)

for i in range(len(r)):
    ax0.scatter(xt[i],0,color=r[i,:],s=100) 

for g,c in zip([gauss1.pdf(np.linspace(-15,15)),gauss2.pdf(np.linspace(-15,15)),gauss3.pdf(np.linspace(-15,15))],['r','g','b']):
    ax0.plot(np.linspace(-15,15),g,color=c,zorder=0)

ax0.set_xlabel('X-axis')
ax0.set_ylabel('Gaussian pdf value')
ax0.legend(['Gaussian 1', 'Gaussian 2', 'Gaussian 3'])
plt.tight_layout()
plt.show()

# M-Step
mc = np.sum(r, axis=0)
pi = mc/len(xt)
mu = np.sum(r*np.vstack((xt, xt, xt)).T, axis=0)/mc
var = []

for c in range(len(pi)):
    var.append(np.sum(np.dot(r[:,c]*(xt[i] - mu[c]).T, r[:,c]*(xt[i] - mu[c])))/mc[c])


class GMM1D:
    """Apply GMM to 1D Data"""
    
    def __init__(self, X, max_iterations):
        """Initialize data and max_iterations"""
        self.X = X
        self.max_iterations = max_iterations
        
    def run(self):
        """Initialize parameters mu, var, pi"""
        self.pi = np.array([1/3, 1/3, 1/3])
        self.mu = np.array([5,8,1])
        self.var = np.array([5,3,1])
        
        r = np.zeros((len(self.X), 3))
        
        for itr in range(self.max_iterations):
    
            gauss1 = norm(loc=self.mu[0], scale=self.var[0])
            gauss2 = norm(loc=self.mu[1], scale=self.var[1])
            gauss3 = norm(loc=self.mu[2], scale=self.var[2])
            
            # E-Step
            for c,g,p in zip(range(3), [gauss1, gauss2, gauss3], self.pi):
                r[:,c] = p*g.pdf(xt[:])

            for i in range(len(r)):
                r[i,:] /= np.sum(r[i,:])

            fig = plt.figure(figsize=(6,5))
            ax0 = fig.add_subplot(111)

            for i in range(len(r)):
                ax0.scatter(xt[i],0,color=r[i,:],s=100) 

            for g,c in zip([gauss1.pdf(np.linspace(-15,15)),gauss2.pdf(np.linspace(-15,15)),gauss3.pdf(np.linspace(-15,15))],['r','g','b']):
                ax0.plot(np.linspace(-15,15),g,color=c,zorder=0)

            plt.tight_layout()
            #plt.show()

            # M-Step
            mc = np.sum(r, axis=0)
            self.pi = mc/len(self.X)
            self.mu = np.sum(r*np.vstack((self.X, self.X, self.X)).T, axis=0)/mc
            self.var = []

            for c in range(len(self.pi)):
                self.var.append(np.sum(np.dot(r[:,c]*(self.X[i] - self.mu[c]).T, r[:,c]*(self.X[i] - self.mu[c])))/mc[c])

gmm = GMM1D(xt, 10)
#gmm.run()

############
# GMM - 2D #
############

from sklearn.datasets import make_blobs
from scipy.stats import multivariate_normal

X,Y = make_blobs(cluster_std=1.5,random_state=20,n_samples=500,centers=3)

X = np.dot(X, np.random.RandomState(0).randn(2,2))
plt.figure(figsize=(6,6))
plt.scatter(X[:, 0], X[:, 1])
plt.show()

class GMM2D:
    """Apply GMM to 2D data"""
    
    def __init__(self, num_clusters, max_iterations):

        """Initialize num_clusters(K) and max_iterations for the model"""

        self.num_clusters = num_clusters
        self.max_iterations = max_iterations

    def run(self, X):

        """Initialize parameters and run E and M step storing log-likelihood value after every iteration"""

        self.pi = np.ones(self.num_clusters)/self.num_clusters
        self.mu = np.random.randint(min(X[:, 0]), max(X[:, 0]), size=(self.num_clusters, len(X[0])))
        self.cov = np.zeros((self.num_clusters, len(X[0]), len(X[0])))

        for n in range(len(self.cov)):
            np.fill_diagonal(self.cov[n], 5)

        # reg_cov is used for numerical stability i.e. to check singularity issues in covariance matrix 
        self.reg_cov = 1e-6*np.identity(len(X[0]))

        x,y = np.meshgrid(np.sort(X[:,0]), np.sort(X[:,1]))
        self.XY = np.array([x.flatten(), y.flatten()]).T
        # Plot the data and the initial model

        fig0 = plt.figure(figsize=(6,6))
        ax0 = fig0.add_subplot(111)
        ax0.scatter(X[:, 0], X[:, 1])
        ax0.set_title("Initial State")

        for m, c in zip(self.mu, self.cov):
            c += self.reg_cov
            multi_normal = multivariate_normal(mean=m, cov=c)
            ax0.contour(np.sort(X[:, 0]), np.sort(X[:, 1]), multi_normal.pdf(self.XY).reshape(len(X), len(X)), colors = 'black', alpha = 0.3)
            ax0.scatter(m[0], m[1], c='grey', zorder=10, s=100)
        
        fig0.savefig('GMM2D Initial State.png')
        plt.show()
        self.log_likelihoods = []

        for iters in range(self.max_iterations):
            # E-Step

            self.ric = np.zeros((len(X), len(self.mu)))

            for pic, muc, covc, r in zip(self.pi, self.mu, self.cov, range(len(self.ric[0]))):
                covc += self.reg_cov
                mn = multivariate_normal(mean=muc, cov=covc)
                self.ric[:, r] = pic*mn.pdf(X)

            for r in range(len(self.ric)):
                self.ric[r, :] = self.ric[r, :] / np.sum(self.ric[r, :])

            # M-step

            self.mc = np.sum(self.ric, axis=0)
            self.pi = self.mc/np.sum(self.mc)
            self.mu = np.dot(self.ric.T, X) / self.mc.reshape(self.num_clusters,1)

            self.cov = []

            for r in range(len(self.pi)):
                covc = 1/self.mc[r] * (np.dot( (self.ric[:, r].reshape(len(X), 1)*(X-self.mu[r]) ).T, X - self.mu[r]) + self.reg_cov)
                self.cov.append(covc)

            self.cov = np.asarray(self.cov)
            self.log_likelihoods.append(np.log(np.sum([self.pi[r]*multivariate_normal(self.mu[r], self.cov[r] + self.reg_cov).pdf(X) for r in range(len(self.pi))])))

            fig1 = plt.figure(figsize=(6,6))
            ax1 = fig1.add_subplot(111)
            ax1.scatter(X[:, 0], X[:, 1])
            ax1.set_title("Iteration " + str(iters))

            for m, c in zip(self.mu, self.cov):
                c += self.reg_cov
                multi_normal = multivariate_normal(mean=m, cov=c)
                ax1.contour(np.sort(X[:, 0]), np.sort(X[:, 1]), multi_normal.pdf(self.XY).reshape(len(X), len(X)), colors = 'black', alpha = 0.3)
                ax1.scatter(m[0], m[1], c='grey', zorder=10, s=100)
            
            fig1.savefig("GMM2D Iter " + str(iters) + ".png")
            plt.show()

        fig2 = plt.figure(figsize=(6,6))
        ax2 = fig2.add_subplot(111)
        ax2.plot(range(0, iters+1, 1), self.log_likelihoods)
        ax2.set_title('Log Likelihood Values')
        fig2.savefig('GMM2D Log Likelihood.png')
        plt.show()

    def predict(self, Y):

        """Predicting cluster for new samples in array Y"""

        predictions = []

        for pic, m, c in zip(self.pi, self.mu, self.cov):
            prob = pic*multivariate_normal(mean=m, cov=c).pdf(Y)
            predictions.append([prob])

        predictions = np.asarray(predictions).reshape(len(Y), self.num_clusters)
        predictions = np.argmax(predictions, axis=1)

        fig2 = plt.figure(figsize=(6,6))
        ax2 = fig2.add_subplot(111)
        ax2.scatter(X[:, 0], X[:, 1], c='c')
        ax2.scatter(Y[:, 0], Y[:, 1], marker='*', c='k', s=150, label = 'New Data')
        ax2.set_title("Predictions on New Data")

        colors = ['r', 'b', 'g']

        for m, c, col, i in zip(self.mu, self.cov, colors, range(len(colors))):
    #         c += reg_cov
            multi_normal = multivariate_normal(mean=m, cov=c)
            ax2.contour(np.sort(X[:, 0]), np.sort(X[:, 1]), multi_normal.pdf(self.XY).reshape(len(X), len(X)), colors = 'black', alpha = 0.3)
            ax2.scatter(m[0], m[1], marker='o', c=col, zorder=10, s=150, label = 'Centroid ' + str(i+1))

        for i in range(len(Y)):
            ax2.scatter(Y[i, 0], Y[i, 1], marker='*', c=colors[predictions[i]], s=150)

        ax2.set_xlabel('X-axis')
        ax2.set_ylabel('Y-axis')
        ax2.legend()
        fig2.savefig('GMM2D Predictions.png')
        plt.show()

        return predictions

y = np.random.randint(-10, 20, size=(12, 2))
gmm2d = GMM2D(num_clusters=3, max_iterations=10)
#gmm2d.run(X)
#gmm2d.predict(y)



from sklearn.mixture import GaussianMixture

X,Y = make_blobs(cluster_std=1.5,random_state=20,n_samples=500,centers=3)
X = np.dot(X, np.random.RandomState(0).randn(2,2))

GMM = GaussianMixture(n_components=3)
GMM.fit(X)
Y = np.random.randint(-10, 20, size=(1, 2))
print(GMM.means_, GMM.predict_proba(Y))

"""Out: 
[[19.88168663 17.47097164] 
[-12.83538784   4.89646199] 
[11.09673732 18.67548935]] 
[[1.91500946e-17 9.30483496e-01 6.95165038e-02]]"""


### Other website
# (https://www.python-course.eu/expectation_maximization_and_gaussian_mixture_models.php)
import matplotlib.pyplot as plt
from matplotlib import style
style.use('fivethirtyeight')
from sklearn.datasets.samples_generator import make_blobs
import numpy as np
from scipy.stats import multivariate_normal
from sklearn.mixture import GaussianMixture

# 0. Create dataset
X,Y = make_blobs(cluster_std=0.5,random_state=20,n_samples=1000,centers=5)

print ('Here')
print (len(X),X.shape,len(Y),Y.shape)
# Stratch dataset to get ellipsoid data
X = np.dot(X,np.random.RandomState(0).randn(2,2))
print (len(X),X.shape)

x,y = np.meshgrid(np.sort(X[:,0]),np.sort(X[:,1]))
XY = np.array([x.flatten(),y.flatten()]).T
print (XY.shape)

GMM = GaussianMixture(n_components=5).fit(X) # Instantiate and fit the model
print('Converged:',GMM.converged_) # Check if the model has converged
means = GMM.means_ 
covariances = GMM.covariances_


# Predict
Y = np.array([[0.5],[0.5]])
prediction = GMM.predict_proba(Y.T)
print(prediction)

# Plot   
fig = plt.figure(figsize=(6,6))
ax0 = fig.add_subplot(111)
ax0.scatter(X[:,0],X[:,1])
ax0.scatter(Y[0,:],Y[1,:],c='orange',zorder=10,s=100)
for m,c in zip(means,covariances):
    multi_normal = multivariate_normal(mean=m,cov=c)
    ax0.contour(np.sort(X[:,0]),np.sort(X[:,1]),multi_normal.pdf(XY).reshape(len(X),len(X)),colors='black',alpha=0.3)
    ax0.scatter(m[0],m[1],c='grey',zorder=10,s=100)
print (len(X),len(multi_normal.pdf(XY)),multi_normal.pdf(XY).shape)
    
plt.show()

####
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn import mixture

n_samples = 300

# generate random sample, two components
np.random.seed(0)

# generate spherical data centered on (20, 20)
shifted_gaussian = np.random.randn(n_samples, 2) + np.array([20, 20])

# generate zero centered stretched Gaussian data
C = np.array([[0., -0.7], [3.5, .7]])
stretched_gaussian = np.dot(np.random.randn(n_samples, 2), C)

# concatenate the two datasets into the final training set
X_train = np.vstack([shifted_gaussian, stretched_gaussian])

# fit a Gaussian Mixture Model with two components
clf = mixture.GaussianMixture(n_components=2, covariance_type='full')
clf.fit(X_train)

# display predicted scores by the model as a contour plot
x = np.linspace(-20., 30.)
y = np.linspace(-20., 40.)
X, Y = np.meshgrid(x, y)
XX = np.array([X.ravel(), Y.ravel()]).T
print (X.shape, Y.shape, XX.shape)
Z = -clf.score_samples(XX)
print (Z.shape)
Z = Z.reshape(X.shape)
print (Z.shape)

CS = plt.contour(X, Y, Z, norm=LogNorm(vmin=1.0, vmax=1000.0),
                 levels=np.logspace(0, 3, 10))
CB = plt.colorbar(CS, shrink=0.8, extend='both')
plt.scatter(X_train[:, 0], X_train[:, 1], .8)

plt.title('Negative log-likelihood predicted by a GMM')
plt.axis('tight')
plt.show()