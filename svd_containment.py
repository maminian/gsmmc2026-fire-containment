import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 14})

# Manuchehr Aminian
# June 11, 2026

####################
# toy model of spread. Can replace this with any callable function which 
# returns X_forecast, an n-by-2 array of the next perimeter.
R0 = 0.5
th = np.arange(0,2*np.pi, 0.1)
x_forecast = lambda tv: (R0+tv)*(np.cos(th) + 0.2*np.cos(5*th))
y_forecast = lambda tv: (R0+tv)*(np.sin(th) + 0.2*np.sin(5*th))
X_forecast = lambda tv: np.vstack([x_forecast(tv), y_forecast(tv)]).T
#############

# containment alg params.
gamma = 2*np.pi +1 # containment speed


dt = 0.1 # right now: time for construction of a "barrier"...?
ts = np.arange(0,2,dt) # for loop and vis; not important

#
num_pts = 3 # number of nearest perimeter points to use for determining directions.
target_dist = 0.4 # ideal distance from local perim (absolute units). idk choice
sigma = 10*gamma*dt # falloff parameter in weighting parallel and perp directions for containment. idk
r = np.zeros((len(ts), 2))
r[0] = [1,0]
dr = np.array([0,1]) # initial direction, or prior direction while in loop. idk if it matters much.

################
#


fig,ax = plt.subplots(constrained_layout=True, figsize=(6,6))

for i_ in range(1,len(ts)):

    #ax.plot(xt,yt, marker='.', c=plt.cm.tab20(i_))
    
    # In general: these are replaced by an array X of points on the fire 
    # perimeter based on any forecasting model.
    t = ts[i_]
    X = X_forecast(t)
    
    # get the nearest point on perimeter and head towards it, a constant number 
    # distance away from the perimeter center
    closest_pts = np.argsort(np.linalg.norm(X-r[i_-1],2,axis=1))

    X_center = np.mean(X, axis=0)
    X_near = X[closest_pts[:num_pts]]
    nearest_bdry = np.mean(X_near, axis=0)
    
    u,s,vh = np.linalg.svd(X_near - nearest_bdry)
    parallel_idx = np.argsort(s) # smallest to largest singular values
    
    #
    parallel_dir = vh[parallel_idx[1],:]
    perp_dir = vh[parallel_idx[0],:]
    
    #
    curr_dist = np.linalg.norm(r[i_-1] - nearest_bdry)

    parallel_vec = np.sign(np.dot(parallel_dir, dr)) * parallel_dir
    # note the minus sign...
    perp_vec = -np.sign(np.dot(perp_dir,r[i_-1] - nearest_bdry)) * perp_dir
    
    
    # magic weighting; want this 1 if we are at the ideal distance; 0 if not close.
    # in theory the perp vec is chosen so it corrects for the direction either way. 
    # outwards if too close; inwards if too far.
    # 
    # I think this weighting is the most subtle part: technically
    # the current version is a 'forward euler' style update; 
    # a step is weighted based on current distance, rather than 
    # what the gap in distances *would be* after a step.
    # TODO: Some kind of rootfinding heuristic would be needed to stabilize.
    alpha = np.exp(-(curr_dist-target_dist)**2/sigma**2)
    
    # identify unit direction of travel based on weighting
    direction = alpha*parallel_vec + np.sqrt(1-alpha**2)*perp_vec
    
    # select next containment point (hence line segment)
    r[i_] = r[i_-1] + dt*gamma*direction
    
    # plot for fun
    if i_==1:
        # label once
        ax.plot(X[:,0],X[:,1], c=plt.cm.tab20(i_-1), alpha=0.2, label='forecasted perimeter')
        ax.scatter(X[closest_pts[:num_pts],0], X[closest_pts[:num_pts],1], c=[plt.cm.tab20(i_)], marker='*', label='reference points')
    else:
        ax.plot(X[:,0],X[:,1], c=plt.cm.tab20(i_-1), alpha=0.2)
        ax.scatter(X[closest_pts[:num_pts],0], X[closest_pts[:num_pts],1], c=[plt.cm.tab20(i_-1)], marker='*')
    

    # hold on to the prior direction for next alg step
    dr = direction
#

ax.scatter(r[:,0], r[:,1], c=plt.cm.tab20(range(len(ts))), marker='s', s=100, label='containment points')
ax.plot(r[:,0], r[:,1], c='k', label='containment')
ax.legend(loc='upper left')

ax.set(aspect='equal', xlabel=r'$x$', ylabel=r'$y$')
fig.show()
fig.savefig('svd_containment_alg.png', bbox_inches='tight')

