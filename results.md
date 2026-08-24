# Ergebnisse

## Hauptkomponentenanalyse (PCA)

Die erste Hauptkomponente (PC1) erklärte 20,58 % und die zweite Hauptkomponente (PC2) 12,26 % der Gesamtvarianz. Zusammen erklärten PC1 und PC2 somit 32,84 % der Varianz im Datensatz. Die Verteilung der erklärten Varianz ist in Abbildung 1 dargestellt.

![Screeplot der PCA](figures/Screeplot.png)

*Abbildung 1: Screeplot der Hauptkomponentenanalyse mit dem Anteil der durch die einzelnen Hauptkomponenten erklärten Varianz.*

In der Darstellung der Proben entlang von PC1 und PC2 zeigte sich keine eindeutige Trennung zwischen den ALL- und AML-Proben. Die beiden Gruppen überlappten über weite Bereiche des PCA-Plots. Einzelne Proben lagen weiter vom Hauptbereich der übrigen Proben entfernt, insgesamt bildeten sich jedoch keine klar voneinander abgegrenzten Cluster entsprechend der Leukämieform. Die Verteilung der ALL- und AML-Proben ist in Abbildung 2 dargestellt.

![PCA der ALL- und AML-Proben](figures/PCA_ALL_AML.png)

*Abbildung 2: PCA der ALL- und AML-Proben anhand der ersten beiden Hauptkomponenten.*

Die Analyse der Loadings zeigte, welche Features besonders stark zu den ersten beiden Hauptkomponenten beitrugen. Zu den stärksten positiven Loadings von PC1 gehörten unter anderem CD47, TGOLN2, CRHBP und TCP1, während BBS9, SMPDL3B, AP3S2 und MPV17 zu den stärksten negativen Loadings zählten. Für PC2 zeigten PRKAR2A, PRKCE, SLC12A5, ZNF44 und NR4A3 starke positive Loadings. Zu den stärksten negativen Loadings gehörten SCN2A, FEV, MAP2K1 und MORF4L2. Die stärksten positiven und negativen Loadings von PC1 sind in Abbildung 3 dargestellt.

![Stärkste Loadings von PC1](figures/PC1_Loadings.png)

*Abbildung 3: Features mit den stärksten positiven und negativen Loadings auf PC1.*