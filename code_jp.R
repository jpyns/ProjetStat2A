rm(list=ls())

donnes<-readRDS("/home/ensai/Bureau/donnee/anticorps_normalisés.rds")


# Transposer le dataframe
data<- as.data.frame(t(donnes))


# Installation des packages nécessaires si ce n'est pas déjà fait
#install.packages("ggplot2")

# Chargement des packages
library(ggplot2)
library(FactoMineR)

# Supposons que 'data' est votre dataframe
# Préparation des données : Standardisation des données
#data_scaled <- scale(data)


head(data)
summary(data)
str(data)


library(corrplot)
library(cowplot)




# Calcul de la matrice de corrélation
correlation_matrix <- cor(data)


# Visualiser la matrice de corrélation
corrplot(correlation_matrix, method = "circle", type = "full", order = "hclust", tl.cex = 0.7)




library(FactoMineR)

# Exécuter l'ACP
res.acp <- PCA(data, graph = TRUE,ncp=4)


library(factoextra)
library(ggplot2)


# Graphique des valeurs propres
fviz_eig(res.acp)


# Graphique des individus
fviz_pca_ind(res.acp)

# Graphique des variables
fviz_pca_var(res.acp)

cols_to_remove<-grep("^Iso",names(data),value=TRUE)
transposed_data_sans_iso<-data[,!names(data)%in%cols_to_remove]

transposed_data_avec_iso<- data[,cols_to_remove]

transposed_data_iso<-data[,cols_to_remove]

res.acp <- PCA(transposed_data_sans_iso, graph = TRUE)


res.acp <- PCA(transposed_data_iso, graph = TRUE)

#faire des tests pour regrouper les isotopes d'apres l'acp precedente trucs de correlation 


#install.packages("cluster", dependencies = TRUE)



library(cluster)

wss <- sapply(1:10, function(k){
  sum(kmeans(transposed_data_sans_iso,
             k, 
             nstart = 10)$withinss)
  })

plot(1:10, wss, type = "b", xlab = "Nombre de Clusters", ylab = "Somme des carrés intra-clusters")




si3<-silhouette(wss[[2]]$cluster, dist = daisy(transposed_data_sans_iso, metric = "euclidean"))

plot(si3)


# Initialiser CPtheta avec une longueur de 10, en ignorant la position 0
CPtheta <- rep(0, 10)
resKmeans <- list()

for (k in 1:10){
  resKmeans[[k]] <- kmeans(transposed_data_sans_iso, k, nstart = 10)
  CPtheta[k] <- resKmeans[[k]]$tot.withinss
}

# Ajuster la commande plot pour ignorer la première valeur de CPtheta (qui est à l'indice 0 et non utilisée)
plot(1:9, CPtheta[-1], type = "b", xlab = "Nombre de Clusters", ylab = "Somme des carrés intra-clusters")


#Pour calculer et visualiser les scores de silhouette en R, 
#qui sont utilisés pour évaluer la qualité des clusters en termes de cohésion et de séparation,
#vous pouvez utiliser la fonction silhouette() disponible dans le package cluster. 
#Le score de silhouette varie de -1 à 1, où une valeur élevée indique que 
#l'objet est bien assorti à son propre cluster et mal assorti aux clusters voisins. Voici comment vous pourriez procéder :



# Supposons que k est le nombre de clusters que vous avez choisi et que vous avez déjà exécuté kmeans
k <- 3  # Exemple, ajustez selon votre cas
resKmeans <- kmeans(transposed_data_sans_iso, k, nstart = 10)



by(transposed_data_sans_iso, resKmeans$cluster,summary)


# Calculer les scores de silhouette
sil <- silhouette(resKmeans$cluster, dist = daisy(transposed_data_sans_iso, metric = "euclidean")


# Visualiser les scores de silhouette
plot(sil, main = paste("Diagramme de silhouette pour k =", k))


by(transposed_data_sans_iso, resKmeans$cluster,summary)

# FAIRE CAH 


library(VarSelLCM)


# Charger le package
library(cluster)

# Exécuter l'analyse hiérarchique
res.CAHward <- agnes(transposed_data_sans_iso, metric = "euclidean", method = "ward")

# Convertir les données en une matrice numérique
data_matrix <- as.matrix(transposed_data_sans_iso)

# Calculer la dissimilarité avec la distance euclidienne
dist_matrix <- dist(data_matrix, method = "euclidean")

# Appliquer la méthode agnes avec "ward" pour le clustering hiérarchique
res.CAHward <- agnes(dist_matrix, method = "ward")

# Afficher les résultats
print(res.CAHward)


res.CAHward <- agnes(transposed_data_sans_iso,metric ="euclidean", method = "ward")
res.CAHmini <- agnes(data,metric ="euclidean", method = "single")
plot(res.CAHward)

# Spécification du nombre de clusters (k)
k <- 2

# Application de l'algorithme K-Means
kmeans_result <- kmeans(transposed_data_sans_iso, centers = k)

# Affichage des résultats
print(kmeans_result)

# Visualisation des clusters
plot(transposed_data_sans_iso, col = kmeans_result$cluster, main = "K-Means Clustering", pch = 16)
points(kmeans_result$centers, col = 1:k, pch = 8, cex = 2)

ggplot(
  data =  kmeans_result$cluster,
  mapping = aes(x = k, y = error)
) +
  geom_line() +
  geom_point() +
  labs(
    x = "k",
    y = "Erreur de validation",
    title = "Evolution de l'erreur selon k"
  
)




library(uwot)
library(umap)
library(reticulate)


# Spécifiez le nombre de dimensions souhaitées (par exemple, n_components = 2 pour une réduction en 2D)
n_components <- 3


res.acp <- PCA(transposed_data_sans_iso, graph = TRUE)

# Exécutez UMAP
resultat_umap <- umap(res.acp n_components = n_components)



##########################

# Préparer les données pour UMAP (utiliser les scores des deux premières composantes, par exemple)
data_for_umap <- res.acp$ind$coord[, 1:n_components] # Assurez-vous que n_components ne dépasse pas le nombre de composantes disponibles

# Exécutez UMAP
resultat_umap <- umap(data_for_umap, n_neighbors = 15, n_components = n_components, min_dist = 0.1)


# Imaginons que vous ayez un vecteur d'étiquettes pour chaque observation
group_labels <- c("Groupe1", "Groupe2", "Groupe1", "groupeff", "fdfdf"), "hjbfejfb")  # Assurez-vous que cela correspond à vos données

# Ensuite, vous pouvez exécuter votre commande de tracé en utilisant 'group_labels' pour colorer
plot(umap_result$layout[,1], umap_result$layout[,2], col = as.factor(group_labels), pch = 16)



# Clustering 

# Utilisation de k-means pour le clustering sur les résultats UMAP
set.seed(123)  # Pour la reproductibilité
km_res <- kmeans(umap_result$layout, centers = 3)  # Ajustez 'centers' au nombre souhaité de clusters



# Installation du package cluster si nécessaire
if(!requireNamespace("cluster", quietly = TRUE)) install.packages("cluster")

# Calcul du score de silhouette
library(cluster)
silhouette_scores <- silhouette(km_res$cluster, dist(umap_result$layout))

# Visualisation du score de silhouette
plot(silhouette_scores, main="Silhouette Plot")




################################
library(Rtsne)

# Spécifiez le nombre de dimensions souhaitées
dims <- 4

# Exécutez t-SNE sur les composantes principales
res.tsne <- Rtsne(res.acp$ind$coord, dims = dims)


# Créez un graphique de dispersion des résultats de t-SNE
plot(res.tsne$Y, col = "blue", pch = 19)


