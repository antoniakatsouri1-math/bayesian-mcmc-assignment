# Bayesian Statistics & MCMC — Εργασία Εξαμήνου

## Περιγραφή

Εργασία εξαμήνου στο μάθημα **Μπεϋζιανή Στατιστική και MCMC**.  
Περιλαμβάνει θεωρητική ανάλυση, υλοποίηση αλγορίθμων MCMC σε R  
και εφαρμογή σε πραγματικά δεδομένα.
---

## Δομή Εργασίας

### Άσκηση 1 — Πολλαπλή Γραμμική Παλινδρόμηση
- Απόδειξη μορφής συνάρτησης πιθανοφάνειας
- Υπολογισμός ύστερων κατανομών (αναλυτικά)
- Εκτίμηση μέσω **JAGS** (Gibbs Sampling)
- Υλοποίηση **Gibbs Sampler** από μηδέν σε R
- Διαγνωστικά MCMC: trace plots, ACF, running quantiles, boxplots

### Άσκηση 2 — Μπεϋζιανός Έλεγχος Υπόθεσης (Poisson)
- Υπολογισμός **Bayes Factor** για H₀: λ=2 vs H₁: λ≠2
- Ύστερη προβλεπτική κατανομή (Αρνητική Διωνυμική)
- Σύγκριση Μπεϋζιανής και κλασικής προσέγγισης

### Άσκηση 3 — Εφαρμογή σε Πραγματικά Δεδομένα (Wish/Kaggle)
- Dataset: *Summer Products Sales Performance* (Kaggle)
- Μοντέλο Binomial με διαφορετικά priors ανά κατηγορία
- Εκτίμηση μέσω **JAGS** με πληροφοριακό και μη-πληροφοριακό prior
- Σύγκριση με κλασική στατιστική (MLE)

### Άσκηση 4 — Jeffreys Prior & Metropolis-Hastings
- Αναλυτική εύρεση Jeffreys prior για Poisson(λ)
- Υλοποίηση **Random Walk Metropolis-Hastings** σε R
- Διαγνωστικά: ergodic mean, ACF, trace plot
- Εκτίμηση 95% Μπεϋζιανού διαστήματος αξιοπιστίας

---

## Τεχνολογίες

- **R** (rjags, R2jags, coda, mvtnorm, MASS, dplyr)
- **JAGS** (Just Another Gibbs Sampler)
- **LaTeX** (Overleaf)

---

## Dataset

[Summer Products Sales Performance — Kaggle](https://www.kaggle.com/datasets/thedevastator/summer-products-sales-performance)

---

## Αρχείο PDF

Το πλήρες κείμενο της εργασίας: `Katsouri-Antonia.pdf`
