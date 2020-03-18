def var_exp(x_obs, y_obs, z_obs):
    vario,dist = nuee_var(x_obs, y_obs, z_obs)
    plt.scatter(dist, vario)

    vario_trie = [x for _, x in sorted(zip(dist, vario))]

    dist_non_trie = dist

    dist.sort()
    print(dist)
    numero_boucle = 0
    valeurs_vario = []
    valeurs_dist = []
    moy_vario = []
    moy_dist = []
    print(len(vario))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in range (math.floor(len(vario)/500)):
            index = 0

            while (index<500):
                valeurs_vario.append(vario_trie[index+numero_boucle*500])
                valeurs_dist.append(dist[index+numero_boucle*500])
                index+=1
            numero_boucle+=1
            moy_vario.append(np.mean(valeurs_vario))
            moy_dist.append(np.mean(valeurs_dist))
            valeurs_vario=[]
            valeurs_dist=[]

        print("Moy vario")
        print(moy_vario)
        print("moy_distance")
        print(moy_dist)

        plt.plot(moy_dist,moy_vario, color='green')

        plt.show()
        return vario,dist_non_trie