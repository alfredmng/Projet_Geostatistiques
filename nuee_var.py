def nuee_var(x_obs, y_obs, z_obs):
    vario = []
    dist = []
    for i in range (len (z_obs)):
        for j in range (len (z_obs)):
            vario.append(0.5 * (z_obs[i] - z_obs[j]) ** 2)
            dist.append(np.sqrt((x_obs[i] - x_obs[j]) ** 2 + (y_obs[i] - y_obs[j]) ** 2))

    plt.scatter(dist,vario)
    plt.show()
    return vario,dist