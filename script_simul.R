a = read.csv2("all_test.csv", sep = ";", header = FALSE, dec = '.')


mutator = c(rep("AID", 10), rep("pmCDA1", 10), rep("rAPOBEC1", 10), rep("TadA*", 10))
a = data.frame(a, mutator = mutator)
co = c(seq(1,11), "mutator")
colnames(a) = co

num = rep(seq(1,11),10)
gg = data.frame(Time = num)
AID = c()
for (i in 31:40) {
  for (j in 1:11)
    AID = c(AID, a[i,j])
}
gg["TadA*"] = AID

ggg = data.frame(Generation = rep(gg$Time,4), 
                  diversity = c(gg$AID, gg$pmCDA1, gg$raPOBEC1, gg$`TadA*`),
                  mutator = c(rep("AID", 110), rep("pmCDA1", 110), rep("rAPOBEC1", 110), rep("TadA*", 110)))

library(ggplot2)

ggplot(ggg) +
  aes(x = Generation, y = diversity, color = mutator) +
  geom_point() +
  geom_smooth()