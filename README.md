# Competitive Contagion of Networks
## Networks
Mathematically, a network or graph $(G)$ can by defined by a set of nodes $(V)$ and the set of links or edges $(E)$ that connect them:
$$G=(V,E).$$

## Competitive contagion
The competitive contagion game involves two players associated with the colors red $(R)$ and blue $(B)$ respectively. At the beginning of the game, each player chooses a node $v \in V$ to influence. Initially, all nodes in the network are in a nuetral state (no color), but those nodes influenced by the players will switch to color $c \in C = \{R,B\}$ with probability:

$$P_c=\frac{\mathbf 1_{c}}{\mathbf 1 _c+\mathbf 1_{C-\{c\}}}.$$

Here, $\mathbf 1_c$ is an indicator function denoting the event where the node has been chosen by the player with color $c$. 

We adopt the generalized adoption function $h(\cdot)$ defined as 
$$h(\alpha_R,\alpha_B)=\frac{\alpha_R}{\alpha_R+\alpha_B}.$$

Notice that this function is competitive since
$$h(\alpha_R,\alpha_B)\leq h(\alpha_A,0)$$

which means that either player benefits from the absence of the other player. 



## The Price of Anarchy (PoA)
In the context of competitive contagion taking place on an arbitrary Network, the Price of Anarchy of the game defined by 
$$\text{PoA}=\frac{\max_{(a_R,a_B)\in \mathcal A_R \times \mathcal A_B }\mathbb E[X_R+X_B|(a_R,a_B)]}{\min_{(\sigma_R,\sigma_B)\in\text{NE}}\mathbb E[X_R+X_B|(\sigma_R,\sigma_B)]}.$$
will be at most $4$. 