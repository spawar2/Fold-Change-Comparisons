# Shrikant Pawar, 07/11/2021

library("readxl")

# Read in the files
MOCK <- read_excel("RNA_data.xlsx", sheet = 1)

KO <- read_excel("RNA_data.xlsx", sheet = 2)

# mocks and KO averages on each replicates for each time-point

mock_Rep1-0 <- mean(MOCK$mock_Rep1_0h, MOCK$mock_Rep2_0h, MOCK$mock_Rep3_0h)

mock_Rep1-2 <- mean(MOCK$mock_Rep1_2h, MOCK$mock_Rep2_2h, MOCK$mock_Rep3_2h)

mock_Rep1-8 <- mean(MOCK$mock_Rep1_8h, MOCK$mock_Rep2_8h, MOCK$mock_Rep3_8h)

mock_Rep1-12 <- mean(MOCK$mock_Rep1_12h, MOCK$mock_Rep2_12h, MOCK$mock_Rep3_12h)


KO_Rep1-0 <- mean(KO$KO_Rep1_0h, KO$KO_Rep2_0h, KO$KO_Rep3_0h)

KO_Rep1-8 <- KO$KO_Rep1_8h, KO$KO_Rep2_8h, KO$KO_Rep3_8h

KO_Rep1-2 <- KO$KO_Rep1_2h, KO$KO_Rep2_2h, KO$KO_Rep3_2h

KO_Rep1-12 <- KO$KO_Rep1_12h, KO$KO_Rep2_12h, KO$KO_Rep3_12h

# Problem 1: How to identify the associated genes and mechanisms with macrophage activations 

# Taking fold changes between averages of KO_Rep1-12/mock_Rep1-12 and KO_Rep1-0/mock_Rep1-0

Fold-change-KO_Rep1-12/mock_Rep1-12 <- KO_Rep1-12/mock_Rep1-12

Fold-change-KO_Rep1-0/mock_Rep1-0 <- KO_Rep1-0/mock_Rep1-0

# Finally fold change between Fold-change-KO_Rep1-12/mock_Rep1-12//Fold-change-KO_Rep1-0/mock_Rep1-0 should give genes significantly upregulated or down regulated with macrophage activations 

macrophage activations <- Fold-change-KO_Rep1-12/mock_Rep1-12/Fold-change-KO_Rep1-0/mock_Rep1-0

# Problem 2: associated genes and mechanisms with macrophage activations in the KO macrophages

# Can be achieved by taking fold changes between KO_Rep1-12/KO_Rep1-0

macrophage activations in the KO <- KO_Rep1-12/KO_Rep1-0

# Problem 3: associated genes and mechanisms with macrophage activations that are changed because of the KO

# This is a difficult problem as there can be several confouunding factors, one strategy can be lokking within timepoint changes by doing following:

fold-KO_Rep1-2/KO_Rep1-0 <- KO_Rep1-2/KO_Rep1-0

fold-KO_Rep1-8/KO_Rep1-2 <- KO_Rep1-8/KO_Rep1-2

fold-KO_Rep1-12/KO_Rep1-8 <- KO_Rep1-12/KO_Rep1-8

# up and down regulated genes at 2, 8 and 12 time-points can be compared to check if there is any pattern/trend of genes changed because of the KO at each of the time-points.




