library(ggplot2)

name = c()
pos = c()
ref = c()
var = c()
read_mut = c()
read_tot = c()
file=list.files(pattern = "A122_10")

for (i in file){
  tmp = read.csv2(i, sep = "\t")
  name = c(name, rep(substring(i, 1,11), nrow(tmp)))
  pos = c(pos, tmp$Position)
  ref = c(ref, tmp$Ref)
  var = c(var, tmp$VarAllele)
  read_mut = c(read_mut, tmp$Reads2)
  read_tot = c(read_tot, tmp$Reads1 + tmp$Reads2)
}

data = data.frame(Sample = name, Position = pos, Ref = ref, Mutation = var,  Occurences = read_mut, Depth = read_tot)
data = data[data$Ref %in% c("A", "T", "C", "G"),]
#write.csv2(data, file = "donneees_ngs_variantcall.csv", sep = ";", dec = ",")

data$Sample = as.factor(data$Sample)
gg = data.frame(data[((data$Ref %in% c('C', 'G')) & data$Mutation %in% c("T", "A") & data$Depth > 20),])

mutators = c("T7", "AID-T7", "pmCDAI-T7", "rAPOBEC1-T7", 
             "TadA*-T7", "CGG", "AID-CGG", "pmCDA1-CGG", "rAPOBEC1-CGG",
             "TadA*-CGG", "evoAPOBEC1-BE4max-CGG", "evo-CDA1-BE4max-CGG",
             "ABE8.20-m-CGG", "evoAPOBEC1-BE4max-T7", "evo-CDA1-BE4max-T7",
             "ABE8.20-m-T7", "221", "pSEVA471")
id = c("100", "101", "102", "103", "104", "109", "110", "111", "112",
       "113", "114", "115", "116", "117", "118", "119", "221", "471")
trans = data.frame(id = id, mut = mutators)

nom = c()
for (i in 1:nrow(gg)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(gg$Sample[i]), 9, 11)])
}
gg[["nom"]] = nom

ggplot(gg) +
  aes(x = nom, y = Occurences/Depth) +
  ggtitle("C:G -> T:A") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gg2 = data
nom = c()
type = c()
for (i in 1:nrow(gg2)) {
  nom = c(nom, trans$mut[trans$id == substring(as.character(gg2$Sample[i]), 9, 11)])
  if (((data$Ref[i] %in% c('C', 'G')) & data$Mutation[i] %in% c("T", "A"))){
    type = c(type, "C:G -> T:A")
  }
  else {
    type = c(type, "A:T -> G:C")
  }
}
gg2[["nom"]] = nom
gg2[["type"]] = type

ggplot(gg2) +
  aes(x = nom, fill = type) +
  ggtitle("total mutations") +
  xlab("Mutator") +
  ylab("mutations") +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

id = levels(as.factor(gg2$nom))
mut = c()
idd = c()
type = c()
for (i in id){
  mut = c(mut, sum(gg2$Occurences[gg2$type == "C:G -> T:A" & gg2$nom == i])/sum(gg2$Depth[gg2$type == "C:G -> T:A" & gg2$nom == i]))
  mut = c(mut, sum(gg2$Occurences[gg2$type == "A:T -> G:C" & gg2$nom == i])/sum(gg2$Depth[gg2$type == "A:T -> G:C" & gg2$nom == i]))
  type = c(type, "C:G -> T:A", "A:T -> G:C")
  idd = c(idd, i, i)
}
gg3 = data.frame(id =idd, mut = mut, type = type)

ggplot(gg3) +
  aes(x = id, y = mut, fill = type) + 
  ggtitle("Mutations of each type") +
  xlab("Mutator") +
  ylab("Proportion") +
  geom_col() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggplot(gg2[gg2$type == "A:T -> G:C",]) +
  aes(x = nom, y = Occurences/Depth) +
  ggtitle("A:T -> G:C") +
  xlab("Mutator") +
  ylab("variant read / total") +
  geom_jitter() +
  geom_boxplot(alpha = 0.5) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))


