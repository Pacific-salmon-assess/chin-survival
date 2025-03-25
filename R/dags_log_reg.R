## Logistic Regression DAGs
# Jun 16, 2024
# Explore conditional dependence of variables influencing Chinook salmon 
# survival

library(dagitty)

dag_simple <-dagitty( 
  "dag {
  L <- D <- A -> S
  A -> F -> S
  A -> L -> S
  D -> F
  D -> S
  S <- Y -> D -> F
  L <- Y -> F
  F -> L 
  }")
adjustmentSets(dag_simple, exposure = "L", outcome = "S")
adjustmentSets(dag_simple, exposure = "F", outcome = "S")
adjustmentSets(dag_simple, exposure = "D", outcome = "S")
adjustmentSets(dag_simple, exposure = "A", outcome = "S")
impliedConditionalIndependencies(dag_simple)

d_model <- bf(D ~ 1 + A)
f_model <- bf(FF ~ 1 + A + Y + D)
l_model <- bf(L ~ 1 + A + Y + D + FF)
s_model <- bf(S ~ 1 + A + Y + D + FF + L)


dag_full <-dagitty( 
  "dag {
  L <- D <- A -> S
  A -> F -> S
  A -> L -> S
  D -> F
  D -> S
  S <- Y -> D -> F
  L <- Y -> F
  F -> L 
  F -> I -> S
  F -> H -> S
  H -> I
  }")
adjustmentSets(dag_full, exposure = "L", outcome = "S")
adjustmentSets(dag_full, exposure = "I", outcome = "S")


# Following 14.4, which parameters should be correlated?
# stock effects on survival, size, capture date, and lipid
# year effects on survival, size, and lipid
# size and lipid effects on survival 

## generative model
online_dag <- dagitty('
dag {
"fork length" [exposure,pos="0.106,0.772"]
"observed survival" [outcome,pos="0.554,1.080"]
"tagging date" [exposure,pos="-0.713,1.122"]
"true survival" [latent,pos="0.374,1.087"]
condition [latent,pos="-0.088,0.591"]
exploitation [exposure,pos="0.307,0.102"]
injury [exposure,pos="0.016,1.481"]
lipid [exposure,pos="0.443,0.574"]
stock [adjusted,pos="-0.734,0.399"]
year [adjusted,pos="-0.230,-0.056"]
"fork length" -> "true survival"
"fork length" -> injury
"tagging date" -> "true survival"
"tagging date" -> condition
"true survival" -> "observed survival"
condition -> "fork length"
condition -> lipid
exploitation -> "true survival"
injury -> "true survival"
lipid -> "true survival"
stock -> "fork length"
stock -> "tagging date"
stock -> "true survival"
stock -> lipid
year -> "true survival"
year -> condition
}
')

## recovery model
online_dag2 <- dagitty('
dag {
"apparent survival" [outcome,pos="0.619,1.400"]
"detection probability" [pos="-0.770,0.090"]
"fork length" [exposure,pos="0.106,0.772"]
"tagging date" [exposure,pos="-0.713,1.122"]
condition [latent,pos="-0.066,0.528"]
exploitation [exposure,pos="0.307,0.102"]
injury [exposure,pos="-0.012,1.599"]
lipid [exposure,pos="0.378,0.688"]
stock [adjusted,pos="-0.860,0.500"]
year [adjusted,pos="-0.244,-0.053"]
"detection probability" -> "apparent survival"
"fork length" -> "apparent survival"
"tagging date" -> "apparent survival"
"tagging date" -> condition
condition -> "fork length"
condition -> lipid
exploitation -> "apparent survival"
injury -> "apparent survival"
lipid -> "apparent survival"
stock -> "apparent survival"
stock -> "tagging date"
stock -> condition
year -> condition
}
'
)


# recover parameters of direct interest; include effects of stock/year to test
# for covariance among intercepts (i.e. stocks with positive condition tend to 
# be captured later in the year and have higher survival); ignore stock/year 
# effects on probability/CYER because including them doesn't impact estimates
# and not of interest causally



## MSF DAG

msf_dag <- dagitty('
dag {
  "fin damage" [pos="-0.032,1.122"]
  "fork length" [exposure,pos="-0.179,0.891"]
  "handling time" [exposure,pos="0.113,0.824"]
  "hook location" [pos="-0.246,1.232"]
  "hook size" [pos="-0.349,1.404"]
  "maturation stage" [pos="0.071,0.069"]
  "net use" [pos="-0.335,1.057"]
  "observed survival" [outcome,pos="0.467,1.071"]
  "scale loss" [pos="0.079,1.039"]
  "tagging date" [exposure,pos="-0.285,0.570"]
  "true survival" [latent,pos="0.279,1.077"]
  SST [exposure,pos="0.242,0.743"]
  condition [latent,pos="-0.088,0.591"]
  exploitation [exposure,pos="0.307,0.102"]
  eye [pos="-0.121,1.437"]
  injury [exposure,pos="0.123,1.258"]
  injury1 [pos="-0.095,1.264"]
  lipid [exposure,pos="0.443,0.574"]
  region [exposure,pos="0.425,0.191"]
  sex [exposure,pos="0.144,0.434"]
  stock [adjusted,pos="-0.579,0.420"]
  year [adjusted,pos="-0.230,-0.056"]
  "fin damage" -> injury
  "fork length" -> "fin damage"
  "fork length" -> "handling time"
  "fork length" -> "hook location"
  "fork length" -> "scale loss"
  "fork length" -> "true survival"
  "handling time" -> "true survival"
  "hook location" -> eye
  "hook location" -> injury1
  "hook size" -> "hook location"
  "hook size" -> eye
  "hook size" -> injury1
  "maturation stage" -> "fin damage"
  "maturation stage" -> "scale loss"
  "maturation stage" -> condition
  "maturation stage" -> sex
  "net use" -> "fin damage"
  "net use" -> "scale loss"
  "scale loss" -> injury
  "tagging date" -> "true survival"
  "tagging date" -> SST
  "tagging date" -> condition
  "true survival" -> "observed survival"
  SST -> "true survival"
  condition -> "fork length"
  condition -> lipid
  exploitation -> "true survival"
  eye -> injury
  injury -> "true survival"
  injury1 -> injury
  lipid -> "true survival"
  region -> "true survival"
  sex -> "true survival"
  stock -> "fork length"
  stock -> "tagging date"
  stock -> "true survival"
  stock -> lipid
  year -> "true survival"
  year -> condition
}
'
)