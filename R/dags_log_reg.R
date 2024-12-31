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
"apparent survival" [outcome,pos="0.665,1.409"]
"detection probability" [pos="-0.770,0.090"]
"fork length" [exposure,pos="0.106,0.772"]
"tagging date" [exposure,pos="-0.713,1.122"]
condition [latent,pos="-0.103,0.591"]
exploitation [exposure,pos="0.307,0.102"]
injury [exposure,pos="-0.012,1.599"]
lipid [exposure,pos="0.378,0.688"]
stock [adjusted,pos="-0.869,0.298"]
year [adjusted,pos="-0.397,-0.024"]
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
stock -> "detection probability"
stock -> "tagging date"
stock -> condition
stock -> exploitation
year -> "apparent survival"
year -> "detection probability"
year -> condition
year -> exploitation
}
')

## recovery model
online_dag2 <- dagitty('
dag {
"apparent survival" [outcome,pos="0.665,1.409"]
"detection probability" [pos="-0.770,0.090"]
"fork length" [exposure,pos="0.106,0.772"]
"tagging date" [exposure,pos="-0.713,1.122"]
condition [latent,pos="-0.103,0.591"]
exploitation [exposure,pos="0.307,0.102"]
injury [exposure,pos="-0.012,1.599"]
lipid [exposure,pos="0.378,0.688"]
stock [adjusted,pos="-0.869,0.298"]
year [adjusted,pos="-0.397,-0.024"]
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
year -> "apparent survival"
}
'
)
