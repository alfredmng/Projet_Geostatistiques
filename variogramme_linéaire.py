def var_lin(x_obs, y_obs, z_obs):
    varexp,dist = var_exp(x_obs, y_obs, z_obs)
    hmax =120 #modifiable
    idx = np.where(np.array(dist) <= hmax)[0]

    #pas besoin si déjà des array
    varexp_array = np.array(varexp)
    dist_array = np.array(dist)

    varexp_array = varexp_array[idx]
    dist_array = dist_array[idx]

    # on force les dimensions en (n,1)
    Y = varexp_array.reshape((-1, 1))
    A = dist_array.reshape((-1, 1))
    Xhat = np.linalg.pinv(A.T @ A) @ A.T @ Y

    var_lin = []
    for i in range (len(dist)):
        var_lin.append(Xhat[0]*dist[i])
    print(var_lin)

    plt.plot(dist, var_lin, color='green')
    plt.show()
