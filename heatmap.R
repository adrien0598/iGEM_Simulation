a = read.csv2("all_test_double.csv", sep = ";", dec = ",", header = FALSE)

library(plotly)

g6 = a[,6]
m = matrix(nrow = 8, ncol = 8, data = g6)
colnames(m) = c("none", "AID", "pmCDA1", "rAPOBEC1", "TadA*", "evoAPOBEC1", "evoCDA1", "ABE.20-m")
row.names(m) = c("none", "AID", "pmCDA1", "rAPOBEC1", "TadA*", "evoAPOBEC1", "evoCDA1", "ABE.20-m")
library(ggplot2)

v1 = rep(c("none", "AID", "pmCDA1", "rAPOBEC1", "TadA*", "evoAPOBEC1", "evoCDA1", "ABE.20-m"), 8)
v2 = c(rep("none", 8), rep("AID", 8), rep("pmCDA1", 8), 
       rep("rAPOBEC1", 8), rep("TadA*", 8), rep("evoAPOBEC1", 8),
       rep("evoCDA1", 8), rep("ABE.20-m", 8))
gg = data.frame(T7 = v1,
                CGG_R12_KIRV = v2)

tmp = c()
for (i in 1:8){
  tmp = c(tmp, as.numeric(m[i,]))
}
gg[["Diversity"]] = tmp
gg <- gg %>% mutate(text = paste0("T7 + ", T7, "\n", "CGG + ", CGG_R12_KIRV, "\n", "Rate: ",round(Mutation_rate,2)))
library(plotly)
g = ggplot(gg) + 
  aes(x = T7, y = CGG_R12_KIRV, fill = Mutation_rate, text = text) +
  scale_x_discrete(limits = c("none", "AID", "pmCDA1", "rAPOBEC1", "TadA*", "evoAPOBEC1", "evoCDA1", "ABE.20-m")) +
  scale_y_discrete(limits = c("none", "AID", "pmCDA1", "rAPOBEC1", "TadA*", "evoAPOBEC1", "evoCDA1", "ABE.20-m")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient2(low = "grey", mid = "orange", high = "red", midpoint = 0.5) +
  geom_tile()

ggplotly(g, tooltip = "text")
install.packages("htmlwidget")
library(htmlwidgets)
saveWidget(pp, file=paste0( getwd(), "/HtmlWidget/ggplotlyHeatmap.html"))
  