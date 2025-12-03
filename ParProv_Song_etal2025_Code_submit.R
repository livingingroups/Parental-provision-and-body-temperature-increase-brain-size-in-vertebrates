#############################
###20250822_by_Zitan#########
#############################
#loading
library(ape)
library(scales)
library(MCMCglmm)
library(phylolm)
library(phytools)
library(geiger)
library(nlme)
library(evomap)
library(fishtree)
#Data input
#ParProv_Song_2025 <- readxl::read_excel( "~/ParProv_Song_etal2025_Data_submit.xlsx", sheet = 1)
rownames(ParProv_Song_2025)<-ParProv_Song_2025$Species
#phylogenic tree input

jawlessTree<-read.tree("~/jawlessTree.tre")
jawlessTree # 16 species

SharkTree<-read.tree("~/SharkTree.tre")
SharkTree # 120 species

fishtree<-fishtree::fishtree_phylogeny()
fishtree<- drop.tip(fishtree, fishtree$tip.label[! fishtree$tip.label %in% ParProv_Song_2025[which( !is.na(ParProv_Song_2025$Residual_BS)),]$Species], root.edge = F, rooted = is.rooted(fishtree))
fishtree # 210 species

AmphibiaTree<-read.tree("~/AmphibiaTree.tre")
AmphibiaTree

CrocoTurtleMCCTree<-read.nexus("~/CrocoTurtleReptileMCCTree.tre") 
CrocoTurtleMCCTree # 221 species

Birdtree<-read.nexus("~/BirdTree2025_1136.tre")
Birdtree #1136 species

MammalTree<-ape::read.tree("~/MammalTree.tre") 
MammalTree # 855 species

AllVerTree<-read.tree("~/AllVerTree.tre") # all vertebrate phylogentice tree
AllVerTree # 2660 species
####

my.prior <- list(R=list(V=1,nu=0.002),
                 G=list(
                   G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100000)))

summary(ParProv_Song_2025$Offspring)
summary(ParProv_Song_2025$BrainMass.BrM..g.)
summary(ParProv_Song_2025$Bodymass.BdM..g.)


##Table 1####

#Lampreys
# Brain size ~ Offspring size 
inv.phylojawlessTree <- inverseA(chronoMPL(jawlessTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025LampreysMCMM_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Lampreys"  )  ,]
nrow(ParProv_Song_2025LampreysMCMM_EA) # 14
ParProv_Song_2025LampreysMCMM_EA$Brain<-NA
ParProv_Song_2025LampreysMCMM_EA$Brain<-scale(log10(ParProv_Song_2025LampreysMCMM_EA$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025LampreysMCMM_EA$Offspring<-scale(log10(ParProv_Song_2025LampreysMCMM_EA$Offspring) , center = T)
ParProv_Song_2025LampreysMCMM_EA$Body<-NA
ParProv_Song_2025LampreysMCMM_EA$Body<-scale(log10(ParProv_Song_2025LampreysMCMM_EA$Bodymass.BdM..g.) , center = T) 

LampreyMCMCALL_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                           random = ~Species , data = ParProv_Song_2025LampreysMCMM_EA, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylojawlessTree$Ainv))
summary(LampreyMCMCALL_log)


VCV<-NA
VCV <- LampreyMCMCALL_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Cartilaginous fishes
# Brain size ~ Offspring size 

ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes" )  ,]
nrow(ParProv_Song_2025SharkMCMC) # 120
ParProv_Song_2025SharkMCMC$Brain<-NA
ParProv_Song_2025SharkMCMC$Brain<-scale(log10(ParProv_Song_2025SharkMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025SharkMCMC$Offspring<-scale(log10(ParProv_Song_2025SharkMCMC$Offspring) , center = T)
ParProv_Song_2025SharkMCMC$Body<-NA
ParProv_Song_2025SharkMCMC$Body<-scale(log10(ParProv_Song_2025SharkMCMC$Bodymass.BdM..g.) , center = T)

inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)

sharkMCMCAll_log<-MCMCglmm(  fixed = Brain~ Body  +Offspring  ,
                         random = ~Species , data = ParProv_Song_2025SharkMCMC, prior = my.prior, 
                         nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCAll_log)

VCV<-NA
VCV <- sharkMCMCAll_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])

mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Ray-finned fishes
# Brain size ~ Offspring size 

inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes" )  ,]
nrow(ParProv_Song_2025fishMCMC) # 210
ParProv_Song_2025fishMCMC$Brain<-NA
ParProv_Song_2025fishMCMC$Brain<-scale(log10(ParProv_Song_2025fishMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025fishMCMC$Offspring<-scale(log10(ParProv_Song_2025fishMCMC$Offspring) , center = T)
ParProv_Song_2025fishMCMC$Body<-NA
ParProv_Song_2025fishMCMC$Body<-scale(log10(ParProv_Song_2025fishMCMC$Bodymass.BdM..g.) , center = T)

fishMCMCAll_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                        random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMCAll_log)


VCV<-NA
VCV <- fishMCMCAll_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Amphibians
# Brain size ~ Offspring size 
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025AmphibainMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC) #  130
ParProv_Song_2025AmphibainMCMC$Brain<-NA
ParProv_Song_2025AmphibainMCMC$Brain<-scale(log10(ParProv_Song_2025AmphibainMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025AmphibainMCMC$Offspring<-scale(log10(ParProv_Song_2025AmphibainMCMC$Offspring) , center = T)
ParProv_Song_2025AmphibainMCMC$Body<-NA
ParProv_Song_2025AmphibainMCMC$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC$Bodymass.BdM..g.) , center = T)

AmphibainMCMCAll_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                             random = ~Species , data = ParProv_Song_2025AmphibainMCMC, prior = my.prior, 
                             nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMCAll_log)

VCV<-NA
VCV <- AmphibainMCMCAll_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Reptiles
# Brain size ~ Offspring size 
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Reptiles" & !is.na(ParProv_Song_2025$Residual_BS)   )  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 190
ParProv_Song_2025ReptileMCMC$Brain<-NA
ParProv_Song_2025ReptileMCMC$Brain<-scale(log10(ParProv_Song_2025ReptileMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025ReptileMCMC$Offspring<-scale(log10(ParProv_Song_2025ReptileMCMC$Offspring) , center = T)
ParProv_Song_2025ReptileMCMC$Body<-NA
ParProv_Song_2025ReptileMCMC$Body<-scale(log10(ParProv_Song_2025ReptileMCMC$Bodymass.BdM..g.) , center = T)

ReptileMCMCAll<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                           random = ~Species , data = ParProv_Song_2025ReptileMCMC, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMCAll)


VCV<-NA
VCV <- ReptileMCMCAll$VCV
lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Birds
# Brain size ~ Offspring size 
inv.phyloBirdTree <- inverseA(chronoMPL(Birdtree) ,nodes="TIPS",scale=TRUE)

ParProvBirdAllMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Birds"    )  ,]
nrow(ParProvBirdAllMCMC) # 1136

ParProvBirdAllMCMC$Brain<-NA
ParProvBirdAllMCMC$Brain<-scale(log10(ParProvBirdAllMCMC$BrainMass.BrM..g.) , center = T)
ParProvBirdAllMCMC$Offspring<-scale(log10(ParProvBirdAllMCMC$Offspring) , center = T)
ParProvBirdAllMCMC$Body<-NA
ParProvBirdAllMCMC$Body<-scale(log10(ParProvBirdAllMCMC$Bodymass.BdM..g.) , center = T)

BirdMCMCAll_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                        random = ~Species , data = ParProvBirdAllMCMC, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree$Ainv))
summary(BirdMCMCAll_log)

VCV<-NA
VCV <- BirdMCMCAll_log$VCV
lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#Mammals
# Brain size ~ Offspring size 
inv.phyloMammTree <- inverseA(chronoMPL(MammalTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025Mamm<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals" ),]
nrow(ParProv_Song_2025Mamm) #  855

ParProv_Song_2025Mamm$Brain<-scale(log10(ParProv_Song_2025Mamm$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025Mamm$Offspring<-scale(log10(ParProv_Song_2025Mamm$Offspring) , center = T)
ParProv_Song_2025Mamm$Body<-NA
ParProv_Song_2025Mamm$Body<-scale(log10(ParProv_Song_2025Mamm$Bodymass.BdM..g.) , center = T)

MammMCMC855_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                        random = ~Species , data = ParProv_Song_2025Mamm, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammTree$Ainv))

summary(MammMCMC855_log)


VCV<-NA
VCV <- MammMCMC855_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])

mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#####


##Table 2####
##all vertebrates

ParProv_Song_2025_All<-ParProv_Song_2025[which(!is.na(ParProv_Song_2025$Tb) ), 
                                         c("Species", "Class", "AdultMass", "BrainMass.BrM..g.", "Offspring", "Tb")]
nrow(ParProv_Song_2025_All) #1059
ParProv_Song_2025_All$AdultMass<-scale(log10(ParProv_Song_2025_All$AdultMass), center = T)
ParProv_Song_2025_All$Brain<-scale(log10(ParProv_Song_2025_All$BrainMass.BrM..g. ), center = T)
ParProv_Song_2025_All$Offspring<-scale(log10( ParProv_Song_2025_All$Offspring), center = T)
ParProv_Song_2025_All$Tb<-scale(ParProv_Song_2025_All$Tb, center = T)
summary(ParProv_Song_2025_All)

summary(phylolm(Brain ~AdultMass +Offspring*Tb, data=ParProv_Song_2025_All, phy = AllVerTree, model = "lambda"    ))

inv.phyloAllVerTree <- inverseA(chronoMPL(AllVerTree) ,nodes="TIPS",scale=TRUE)

AllVerMCMC1059_log<-MCMCglmm(  fixed = Brain~  AdultMass +Offspring*Tb    ,
                               random = ~Species , data = ParProv_Song_2025_All, prior = my.prior, 
                               nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAllVerTree$Ainv))

summary(AllVerMCMC1059_log)

VCV<-NA
VCV <- AllVerMCMC1060_log$VCV

lambda_post <- VCV[, "Species"] / (VCV[, "Species"] + VCV[, "units"])
mean_lambda <- mean(lambda_post)
CI_lambda <- quantile(lambda_post, c(0.025, 0.975))
paste(round(mean_lambda, 3), " (", round(CI_lambda[1],3), ", ", round(CI_lambda[2],3), ")", sep ="")

#####


##Table S3####
#Cartilaginous fishes
inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025SharkMCMC) # 120

sharkMCMCOff_Care_log<-MCMCglmm(  fixed = log10(Offspring)~ log10(AdultMass)  +EggCare5Cat  ,
                                  random = ~Species , data = ParProv_Song_2025SharkMCMC, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCOff_Care_log)



sharkMCMCOff_Care_LB_Prov_log<-MCMCglmm(  fixed = log10(Offspring)~ log10(AdultMass)  +EggCare5Cat  ,
                                          random = ~Species , data = ParProv_Song_2025SharkMCMC[which(ParProv_Song_2025SharkMCMC$EggCare5Cat%in% c("3Bear", "4Pre-hat. prov.")),], prior = my.prior, 
                                          nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCOff_Care_LB_Prov_log)

#Ray-finned fishes
inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes" & !is.na(ParProv_Song_2025$EggCare5Cat) )  ,]
nrow(ParProv_Song_2025fishMCMC) # 196

fishMCMC_Off_Care_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                  random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMC_Off_Care_log)


fishMCMC_Off_Care_EG_LB_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                        random = ~Species , data = ParProv_Song_2025fishMCMC[which(ParProv_Song_2025fishMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),], prior = my.prior, 
                                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMC_Off_Care_EG_LB_log)


#Amphibians
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025AmphibainMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC) #  119

AmphibainMCMC_Off_care_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                       random = ~Species , data = ParProv_Song_2025AmphibainMCMC, prior = my.prior, 
                                       nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_Off_care_log)


AmphibainMCMC_Off_care_EG_LB_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                             random = ~Species , data = ParProv_Song_2025AmphibainMCMC[which(ParProv_Song_2025AmphibainMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),], prior = my.prior, 
                                             nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_Off_care_EG_LB_log)

#Reptiles
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)
# Offspring~Egg care
ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Reptiles"    & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 117

ReptileMCMC_Off_Care_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                     random = ~Species , data = ParProv_Song_2025ReptileMCMC  , prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))

summary(ReptileMCMC_Off_Care_log)

ReptileMCMC_Off_Care_EG_LB_log<-MCMCglmm(  fixed = log10(Offspring) ~ log10(AdultMass)  +EggCare5Cat  ,
                                           random = ~Species , data = ParProv_Song_2025ReptileMCMC[which(ParProv_Song_2025ReptileMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),]  , prior = my.prior, 
                                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))

summary(ReptileMCMC_Off_Care_EG_LB_log)

#####


##Table S4####
#Cartilaginous fishes
# Brain size ~ Offspring size 
inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)

# Abandon
ParProv_Song_2025SharkMCMM_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes" & ParProv_Song_2025$EggCare5Cat =="1Abandon" )  ,]
nrow(ParProv_Song_2025SharkMCMM_EA) # 26
ParProv_Song_2025SharkMCMM_EA$Brain<-NA
ParProv_Song_2025SharkMCMM_EA$Brain <-scale(log10(ParProv_Song_2025SharkMCMM_EA$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025SharkMCMM_EA$Offspring<-scale(log10(ParProv_Song_2025SharkMCMM_EA$Offspring) , center = T)
ParProv_Song_2025SharkMCMM_EA$Body<-NA
ParProv_Song_2025SharkMCMM_EA$Body<-scale(log10(ParProv_Song_2025SharkMCMM_EA$Bodymass.BdM..g.) , center = T) 

sharkMCMC_EA_scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                   random = ~Species , data = ParProv_Song_2025SharkMCMM_EA, prior = my.prior, 
                                   nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_EA_scale_log)

#Bear
ParProv_Song_2025SharkMCMM_LB<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"   & ParProv_Song_2025$EggCare5Cat =="3Bear" )  ,]
nrow(ParProv_Song_2025SharkMCMM_LB) # 36
ParProv_Song_2025SharkMCMM_LB$Brain<-NA
ParProv_Song_2025SharkMCMM_LB$Brain<-scale(log10(ParProv_Song_2025SharkMCMM_LB$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025SharkMCMM_LB$Offspring<-scale(log10(ParProv_Song_2025SharkMCMM_LB$Offspring) , center = T)
ParProv_Song_2025SharkMCMM_LB$Body<-NA
ParProv_Song_2025SharkMCMM_LB$Body<-scale(log10(ParProv_Song_2025SharkMCMM_LB$Bodymass.BdM..g.) , center = T) 

sharkMCMC_Bear_scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                     random = ~Species , data = ParProv_Song_2025SharkMCMM_LB, prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_Bear_scale_log)

#Pre-hatching provision
ParProv_Song_2025SharkMCMM_PreH<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & ParProv_Song_2025$EggCare5Cat =="4Pre-hat. prov." )  ,]
nrow(ParProv_Song_2025SharkMCMM_PreH) # 58
ParProv_Song_2025SharkMCMM_PreH$Brain<-NA
ParProv_Song_2025SharkMCMM_PreH$Brain<-scale(log10(ParProv_Song_2025SharkMCMM_PreH$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025SharkMCMM_PreH$Offspring<-scale(log10(ParProv_Song_2025SharkMCMM_PreH$Offspring) , center = T)
ParProv_Song_2025SharkMCMM_PreH$Body<-NA
ParProv_Song_2025SharkMCMM_PreH$Body<-scale(log10(ParProv_Song_2025SharkMCMM_PreH$Bodymass.BdM..g.) , center = T) 

sharkMCMCAllPre_H_scale_log<-MCMCglmm(  fixed = Brain ~Body  +Offspring  ,
                                        random = ~Species , data = ParProv_Song_2025SharkMCMM_PreH, prior = my.prior, 
                                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMCAllPre_H_scale_log)


#Ray-finned fishe
# Brain size ~ Offspring size 
inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025fishMCMC_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & ParProv_Song_2025$EggCare5Cat =="1Abandon"  )  ,]
nrow(ParProv_Song_2025fishMCMC_EA) # 107
ParProv_Song_2025fishMCMC_EA$Brain<-NA
ParProv_Song_2025fishMCMC_EA$Brain<-scale(log10(ParProv_Song_2025fishMCMC_EA$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025fishMCMC_EA$Offspring<-scale(log10(ParProv_Song_2025fishMCMC_EA$Offspring) , center = T)
ParProv_Song_2025fishMCMC_EA$Body<-NA
ParProv_Song_2025fishMCMC_EA$Body<-scale(log10(ParProv_Song_2025fishMCMC_EA$Bodymass.BdM..g.) , center = T)

fishMCMC_EA_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                  random = ~Species , data = ParProv_Song_2025fishMCMC_EA, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_EA_Scale_log)

#Guard
ParProv_Song_2025fishMCMC_EG<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & ParProv_Song_2025$EggCare5Cat =="2Guard" )  ,]
nrow(ParProv_Song_2025fishMCMC_EG) # 42
ParProv_Song_2025fishMCMC_EG$Brain<-NA
ParProv_Song_2025fishMCMC_EG$Brain<-scale(log10(ParProv_Song_2025fishMCMC_EG$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025fishMCMC_EG$Offspring<-scale(log10(ParProv_Song_2025fishMCMC_EG$Offspring) , center = T)
ParProv_Song_2025fishMCMC_EG$Body<-NA
ParProv_Song_2025fishMCMC_EG$Body<-scale(log10(ParProv_Song_2025fishMCMC_EG$Bodymass.BdM..g.) , center = T)

fishMCMC_EG_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                  random = ~Species , data = ParProv_Song_2025fishMCMC_EG, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_EG_Scale_log)


#Bear
ParProv_Song_2025fishMCMC_LB<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"   & ParProv_Song_2025$EggCare5Cat =="3Bear"  )  ,]
nrow(ParProv_Song_2025fishMCMC_LB) # 47
ParProv_Song_2025fishMCMC_LB$Brain<-NA
ParProv_Song_2025fishMCMC_LB$Brain<-scale(log10(ParProv_Song_2025fishMCMC_LB$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025fishMCMC_LB$Offspring<-scale(log10(ParProv_Song_2025fishMCMC_LB$Offspring) , center = T)
ParProv_Song_2025fishMCMC_LB$Body<-NA
ParProv_Song_2025fishMCMC_LB$Body<-scale(log10(ParProv_Song_2025fishMCMC_LB$Bodymass.BdM..g.) , center = T)

fishMCMC_LB_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                  random = ~Species , data = ParProv_Song_2025fishMCMC_LB, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_LB_Scale_log)


#Amphibians
# Brain size ~ Offspring size
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025AmphibainMCMC_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"   & ParProv_Song_2025$EggCare5Cat == "1Abandon"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC_EA) #  90
ParProv_Song_2025AmphibainMCMC_EA$Brain<-NA
ParProv_Song_2025AmphibainMCMC_EA$Brain<-scale(log10(ParProv_Song_2025AmphibainMCMC_EA$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025AmphibainMCMC_EA$Offspring<-scale(log10(ParProv_Song_2025AmphibainMCMC_EA$Offspring) , center = T)
ParProv_Song_2025AmphibainMCMC_EA$Body<-NA
ParProv_Song_2025AmphibainMCMC_EA$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC_EA$Bodymass.BdM..g.) , center = T)

AmphibainMCMC_EA_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                       random = ~Species , data = ParProv_Song_2025AmphibainMCMC_EA, prior = my.prior, 
                                       nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_EA_Scale_log)

#Guard
ParProv_Song_2025AmphibainMCMC_EG<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  & ParProv_Song_2025$EggCare5Cat == "2Guard"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC_EG) #  25
ParProv_Song_2025AmphibainMCMC_EG$Brain<-NA
ParProv_Song_2025AmphibainMCMC_EG$Brain<-scale(log10(ParProv_Song_2025AmphibainMCMC_EG$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025AmphibainMCMC_EG$Offspring<-scale(log10(ParProv_Song_2025AmphibainMCMC_EG$Offspring) , center = T)
ParProv_Song_2025AmphibainMCMC_EG$Body<-NA
ParProv_Song_2025AmphibainMCMC_EG$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC_EG$Bodymass.BdM..g.) , center = T)

AmphibainMCMC_EG_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring ,
                                       random = ~Species , data = ParProv_Song_2025AmphibainMCMC_EG, prior = my.prior, 
                                       nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_EG_Scale_log)


#Reptiles
#Brain size ~ Offspring size
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025ReptileMCMC_EA<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"    & ParProv_Song_2025$EggCare5Cat=="1Abandon" )  ,]
nrow(ParProv_Song_2025ReptileMCMC_EA) # 66
ParProv_Song_2025ReptileMCMC_EA$Brain<-NA
ParProv_Song_2025ReptileMCMC_EA$Brain<-scale(log10(ParProv_Song_2025ReptileMCMC_EA$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025ReptileMCMC_EA$Offspring<-scale(log10(ParProv_Song_2025ReptileMCMC_EA$Offspring) , center = T)
ParProv_Song_2025ReptileMCMC_EA$Body<-NA
ParProv_Song_2025ReptileMCMC_EA$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_EA$Bodymass.BdM..g.) , center = T)


ReptileMCMC_EA_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                     random = ~Species , data = ParProv_Song_2025ReptileMCMC_EA, prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_EA_Scale_log)


#2Guard
ParProv_Song_2025ReptileMCMC_EG<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"    & ParProv_Song_2025$EggCare5Cat =="2Guard"   )  ,]
nrow(ParProv_Song_2025ReptileMCMC_EG) # 21
ParProv_Song_2025ReptileMCMC_EG$Brain<-NA
ParProv_Song_2025ReptileMCMC_EG$Brain<-scale(log10(ParProv_Song_2025ReptileMCMC_EG$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025ReptileMCMC_EG$Offspring<-scale(log10(ParProv_Song_2025ReptileMCMC_EG$Offspring) , center = T)
ParProv_Song_2025ReptileMCMC_EG$Body<-NA
ParProv_Song_2025ReptileMCMC_EG$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_EG$Bodymass.BdM..g.) , center = T)

ReptileMCMC_EG_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                     random = ~Species , data = ParProv_Song_2025ReptileMCMC_EG, prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_EG_Scale_log)


#Bear
ParProv_Song_2025ReptileMCMC_LB<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"   &  ParProv_Song_2025$EggCare5Cat=="3Bear"  )  ,]
nrow(ParProv_Song_2025ReptileMCMC_LB) # 34
ParProv_Song_2025ReptileMCMC_LB$Brain<-NA
ParProv_Song_2025ReptileMCMC_LB$Brain<-scale(log10(ParProv_Song_2025ReptileMCMC_LB$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025ReptileMCMC_LB$Offspring<-scale(log10(ParProv_Song_2025ReptileMCMC_LB$Offspring) , center = T)
ParProv_Song_2025ReptileMCMC_LB$Body<-NA
ParProv_Song_2025ReptileMCMC_LB$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_LB$Bodymass.BdM..g.) , center = T)

ReptileMCMC_LB_Scale_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                     random = ~Species , data = ParProv_Song_2025ReptileMCMC_LB, prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_LB_Scale_log)


#Birds
# Brain size ~ Offspring size
# all incubation
inv.phyloBirdTree <- inverseA(chronoMPL(Birdtree) ,nodes="TIPS",scale=TRUE)

ParProvBirdAllMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Birds"    )  ,]
nrow(ParProvBirdAllMCMC) # 1136
ParProvBirdAllMCMC$Brain<-NA
ParProvBirdAllMCMC$Brain<-scale(log10(ParProvBirdAllMCMC$BrainMass.BrM..g.) , center = T)
ParProvBirdAllMCMC$Offspring<-scale(log10(ParProvBirdAllMCMC$Offspring) , center = T)
ParProvBirdAllMCMC$Body<-NA
ParProvBirdAllMCMC$Body<-scale(log10(ParProvBirdAllMCMC$Bodymass.BdM..g.) , center = T)

BirdMCMC_Pre_Hatch_incubte_log<-MCMCglmm(  fixed = Brain~Body  +Offspring  ,
                                           random = ~Species , data = ParProvBirdAllMCMC[which( ParProvBirdAllMCMC$EggCare5Cat!="1Abandon"),], prior = my.prior, 
                                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree$Ainv))
summary(BirdMCMC_Pre_Hatch_incubte_log)

#####


##Estimate adult brain####
#Adult brain size is estimated by applying a model that is fitted using raw brain and body mass data inculding juveniles. The estimation for adults is made by predicting brain mass based on adult body mass and subsequently adding the residual of the brain mass from the model to retain species-specific trends in brain size relative to body mass.

ParProv_Song_2025$AdultBrain<-NA
#Hagfishes
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Hagfishes"),]$AdultMass_note)
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Hagfishes"),]$AdultBrain<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Hagfishes"),]$BrainMass.BrM..g.
#Lampreys
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultMass_note)
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.) ~ log10(Bodymass.BdM..g.), ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),])
summary(a)  
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$Residual_BS<-a$residuals

for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultBrain[i]<- 10^(predict(a,newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"),]$Residual_BS[i])
  }
}

#Cartilaginous fishes
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultMass_note)
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Cartilaginous fishes" ),] )
summary(a)

for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultBrain[i]<- 10^(predict(a,newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"),]$Residual_BS[i])
  }
}

#Lungfishes
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lungfishes" ),]$AdultMass_note)
a<-NA
b<-ParProv_Song_2025[which( ParProv_Song_2025$Class %in% c("Ray-finned fishes", "Lungfishes") ),]
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class %in% c("Ray-finned fishes", "Lungfishes") ),] )
summary(a)
b$Residual_BS<-a$residuals

for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lungfishes"),])){
  
  ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lungfishes"),]$AdultBrain[i]<- 10^(predict(a,newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lungfishes"),]$AdultMass[i])) + b[which(b$Class=="Lungfishes"),]$Residual_BS[i])
  
}

#Coelacanths
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Coelacanths"  ),]$AdultMass_note)
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Coelacanths"  ),]$AdultBrain<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Coelacanths"  ),]$BrainMass.BrM..g.

#Ray-finned fishes
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),]$AdultMass_note)
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),]$AdultBrain<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),]$BrainMass.BrM..g.

#Amphibians
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"  ),]$AdultMass_note)
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians" ),] )
summary(a)

for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$AdultBrain[i]<- 10^(predict(a,newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"),]$Residual_BS[i])
  }
}

#Reptiles
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"  ),]$AdultMass_note)
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles" & !is.na(ParProv_Song_2025$Offspring)),] )
summary(a)
for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$AdultBrain[i]<- 10^(predict(a, newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"),]$Residual_BS[i])
  }
}


#Birds
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"  ),]$AdultMass_note)
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Birds" & !is.na(ParProv_Song_2025$Offspring) ),] )
summary(a)
for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$AdultBrain[i]<- 10^(predict(a, newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$Residual_BS[i])
  }
}

#Mammals
table(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"  ),]$AdultMass_note)
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Mammals" & !is.na(ParProv_Song_2025$Offspring) ),] )
summary(a)
for(i in 1: nrow(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),])){
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$AdultMass_note[i] == "brain body mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$AdultBrain[i]<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$BrainMass.BrM..g.[i]
  }
  if(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$AdultMass_note[i] == "adult mass"){
    ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$AdultBrain[i]<- 10^(predict(a, newdata = data.frame(Bodymass.BdM..g. = ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$AdultMass[i])) + ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$Residual_BS[i])
  }
}

#####


##Table S5####
Lamprey_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"  ),], phy=jawlessTree, model = "lambda")

shark_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"  ),], phy=SharkTree, model = "lambda")

fish_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),], phy=fishtree, model = "lambda")

Amphi_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"  ),], phy=AmphibiaTree, model = "lambda")

Reptil_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"  ),], phy=CrocoTurtleMCCTree, model = "lambda")

Bird_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"  ),], phy=Birdtree, model = "lambda")

Mammal_Brain_Off<-phylolm(log10(AdultBrain)~log10(Offspring),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"  ),], phy=MammalTree, model = "lambda")

#####


##Table S6####
Lamprey_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Lampreys"  ),], phy=jawlessTree, model = "lambda")

shark_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Cartilaginous fishes"  ),], phy=SharkTree, model = "lambda")

fish_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Ray-finned fishes"  ),], phy=fishtree, model = "lambda")

Amphi_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Amphibians"  ),], phy=AmphibiaTree, model = "lambda")

Reptil_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Reptiles"  ),], phy=CrocoTurtleMCCTree, model = "lambda")

Bird_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"  ),], phy=Birdtree, model = "lambda")

Mammal_brain_body<-phylolm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),  data=ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"  ),], phy=MammalTree, model = "lambda")

#####


##Table S7####

# shark
ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025SharkMCMC) # 118
ParProv_Song_2025SharkMCMC$Brain<-NA
ParProv_Song_2025SharkMCMC$Brain<-scale(log10(ParProv_Song_2025SharkMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025SharkMCMC$Offspring<-scale(log10(ParProv_Song_2025SharkMCMC$Offspring) , center = T)
ParProv_Song_2025SharkMCMC$Body<-NA
ParProv_Song_2025SharkMCMC$Body<-scale(log10(ParProv_Song_2025SharkMCMC$Bodymass.BdM..g.) , center = T)

ParProv_Song_2025SharkMCMC$AdultMass<-scale(log10(ParProv_Song_2025SharkMCMC$AdultMass) , center = T)
ParProv_Song_2025SharkMCMC$Tb<-scale(ParProv_Song_2025SharkMCMC$Tb , center = T)


inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)


sharkMCMC_BS_Tb_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring+ Tb  ,
                                random = ~Species , data = ParProv_Song_2025SharkMCMC, prior = my.prior, 
                                nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_BS_Tb_log)


## fish 
ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025fishMCMC) # 159

ParProv_Song_2025fishMCMC$Brain<-NA
ParProv_Song_2025fishMCMC$Brain<-scale(log10(ParProv_Song_2025fishMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025fishMCMC$Offspring<-scale(log10(ParProv_Song_2025fishMCMC$Offspring) , center = T)

ParProv_Song_2025fishMCMC$Body<-NA
ParProv_Song_2025fishMCMC$Body<-scale(log10(ParProv_Song_2025fishMCMC$Bodymass.BdM..g.) , center = T)

ParProv_Song_2025fishMCMC$AdultMass<-scale(log10(ParProv_Song_2025fishMCMC$AdultMass) , center = T)
ParProv_Song_2025fishMCMC$Tb<-scale(ParProv_Song_2025fishMCMC$Tb , center = T)

inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)



##brain
fishMCMCBS_Tb_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring+ Tb  , 
                              random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMCBS_Tb_log)




## Amphibias 
ParProv_Song_2025AmphibiansMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Amphibians"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025AmphibiansMCMC) # 41
ParProv_Song_2025AmphibiansMCMC$Brain<-NA
ParProv_Song_2025AmphibiansMCMC$Brain<-scale(log10(ParProv_Song_2025AmphibiansMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025AmphibiansMCMC$Offspring<-scale(log10(ParProv_Song_2025AmphibiansMCMC$Offspring) , center = T)

ParProv_Song_2025AmphibiansMCMC$Body<-NA
ParProv_Song_2025AmphibiansMCMC$Body<-scale(log10(ParProv_Song_2025AmphibiansMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025AmphibiansMCMC$AdultMass<-scale(log10(ParProv_Song_2025AmphibiansMCMC$AdultMass) , center = T)

ParProv_Song_2025AmphibiansMCMC$Tb<-scale(ParProv_Song_2025AmphibiansMCMC$Tb , center = T)

inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)


   
AmphibianMCMCBrain_Temp_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring +Tb  ,
                                        random = ~Species , data = ParProv_Song_2025AmphibiansMCMC, prior = my.prior, 
                                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))
summary(AmphibianMCMCBrain_Temp_log)




## Reptile 
ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Reptiles"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 112

ParProv_Song_2025ReptileMCMC$Brain<-NA
ParProv_Song_2025ReptileMCMC$Brain<-scale(log10(ParProv_Song_2025ReptileMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025ReptileMCMC$Offspring<-scale(log10(ParProv_Song_2025ReptileMCMC$Offspring) , center = T)

ParProv_Song_2025ReptileMCMC$Body<-NA
ParProv_Song_2025ReptileMCMC$Body<-scale(log10(ParProv_Song_2025ReptileMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025ReptileMCMC$AdultMass<-scale(log10(ParProv_Song_2025ReptileMCMC$AdultMass) , center = T)

ParProv_Song_2025ReptileMCMC$Tb<-scale(ParProv_Song_2025ReptileMCMC$Tb , center = T)

inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)


ReptileMCMCBrain_Temp_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring +Tb  ,
                                      random = ~Species , data = ParProv_Song_2025ReptileMCMC, prior = my.prior, 
                                      nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMCBrain_Temp_log)


## Birds 
ParProv_Song_2025BirdMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Birds"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025BirdMCMC) # 348

ParProv_Song_2025BirdMCMC$Brain<-NA
ParProv_Song_2025BirdMCMC$Brain<-scale(log10(ParProv_Song_2025BirdMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025BirdMCMC$Offspring<-scale(log10(ParProv_Song_2025BirdMCMC$Offspring) , center = T)

ParProv_Song_2025BirdMCMC$Body<-NA
ParProv_Song_2025BirdMCMC$Body<-scale(log10(ParProv_Song_2025BirdMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025BirdMCMC$AdultMass<-scale(log10(ParProv_Song_2025BirdMCMC$AdultMass) , center = T)

ParProv_Song_2025BirdMCMC$Tb<-scale(ParProv_Song_2025BirdMCMC$Tb , center = T)

BirdTree4<-drop.tip(Birdtree, Birdtree$tip.label[! Birdtree$tip.label %in% ParProv_Song_2025BirdMCMC$Species],
                    root.edge = F, rooted = is.rooted(Birdtree))
BirdTree4
inv.phyloBirdTree4 <- inverseA(chronoMPL(BirdTree4) ,nodes="TIPS",scale=TRUE)


BirdMCMCBrain_Temp_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring +Tb  ,
                                   random = ~Species , data = ParProv_Song_2025BirdMCMC, prior = my.prior, 
                                   nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree4$Ainv))
summary(BirdMCMCBrain_Temp_log)


## Mammals 
ParProv_Song_2025MammalMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Mammals"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025MammalMCMC) # 274

ParProv_Song_2025MammalMCMC$Brain<-NA
ParProv_Song_2025MammalMCMC$Brain<-scale(log10(ParProv_Song_2025MammalMCMC$BrainMass.BrM..g.) , center = T)
ParProv_Song_2025MammalMCMC$Offspring<-scale(log10(ParProv_Song_2025MammalMCMC$Offspring) , center = T)

ParProv_Song_2025MammalMCMC$Body<-NA
ParProv_Song_2025MammalMCMC$Body<-scale(log10(ParProv_Song_2025MammalMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025MammalMCMC$AdultMass<-scale(log10(ParProv_Song_2025MammalMCMC$AdultMass) , center = T)

ParProv_Song_2025MammalMCMC$Tb<-scale(ParProv_Song_2025MammalMCMC$Tb , center = T)


MammalTree1<-drop.tip(MammalTree, MammalTree$tip.label[! MammalTree$tip.label %in% ParProv_Song_2025MammalMCMC$Species],
                      root.edge = F, rooted = is.rooted(MammalTree))
inv.phyloMammTree1 <- inverseA(chronoMPL(MammalTree1) ,nodes="TIPS",scale=TRUE)



MammalMCMCBrain_Temp_log<-MCMCglmm(  fixed = Brain ~  Body +Offspring +Tb  ,
                                     random = ~Species , data = ParProv_Song_2025MammalMCMC, prior = my.prior, 
                                     nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammTree1$Ainv))
summary(MammalMCMCBrain_Temp_log)


#####



##Table S8-S10, BMR in Mammals & Birds####
ParProv_Song_2025$Residual_bmr<-NA

ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$Residual_bmr<-NA
a<-NA
a<-lm(log10(BMR) ~ log10(BMR_mass), ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),])
summary(a)  
for(i in 1 :length(a$residuals)){
  ParProv_Song_2025[which(ParProv_Song_2025$Species ==names(a$residuals[i]) ),]$Residual_bmr<-as.numeric(a$residuals[i])
  
}

 

ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$Residual_bmr<-NA
a<-NA
a<-lm(log10(BMR) ~ log10(BMR_mass), ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),])
summary(a)  
for(i in 1 :length(a$residuals)){
  ParProv_Song_2025[which(ParProv_Song_2025$Species ==names(a$residuals[i]) ),]$Residual_bmr<-as.numeric(a$residuals[i])
  
}

ParProv_Song_2025$BMQ_scale<-NA
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$BMQ_scale<-(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$BMR*86.4)/(299*(1.05*(ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$BMR_mass/1000)*0.95)^0.72*10^(-0.0057*ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$Temp))

#Mammal

ParProv_Song_2025PathMammal_BMR<-ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Mammals")
                                                         & !is.na(ParProv_Song_2025$Tb) 
                                                         & !is.na(ParProv_Song_2025$Residual_bmr) 
                                                         & !is.na(ParProv_Song_2025$Temp) ),
                                                   c( "Species",   "Residual_BS",
                                                      "Residual_OS",
                                                      "Bodymass.BdM..g.",
                                                      "Residual_bmr", 
                                                      "Tb","Temp" ,"BMQ_scale"
                                                   )]

nrow(ParProv_Song_2025PathMammal_BMR)
summary(ParProv_Song_2025PathMammal_BMR)

ParProv_Song_2025PathMammal_BMR$Brain<-scale((ParProv_Song_2025PathMammal_BMR$Residual_BS), center = T)

ParProv_Song_2025PathMammal_BMR$Newborn <-scale((ParProv_Song_2025PathMammal_BMR$Residual_OS), center = T)
ParProv_Song_2025PathMammal_BMR$Body<-scale(log10(ParProv_Song_2025PathMammal_BMR$Bodymass.BdM..g.), center = T)

ParProv_Song_2025PathMammal_BMR$BMR <-scale((ParProv_Song_2025PathMammal_BMR$Residual_bmr), center = T)
ParProv_Song_2025PathMammal_BMR$BMQ <-scale((ParProv_Song_2025PathMammal_BMR$BMQ_scale), center = T)

ParProv_Song_2025PathMammal_BMR$Tb<-scale((ParProv_Song_2025PathMammal_BMR$Tb), center = T)
ParProv_Song_2025PathMammal_BMR$Ta<-scale((ParProv_Song_2025PathMammal_BMR$Temp), center = T)

summary(ParProv_Song_2025PathMammal_BMR)


MammalTree1 <- drop.tip(MammalTree , MammalTree$tip.label[which(MammalTree$tip.label %ni% rownames(ParProv_Song_2025PathMammal_BMR))])
MammalTree1
inv.phyloMammalTree1 <- inverseA(chronoMPL(MammalTree1) ,nodes="TIPS",scale=TRUE)

MammalMCMC_BMR<-MCMCglmm(  fixed = Brain~Body  +Newborn + BMR + Tb +Ta,
                           random = ~Species , data = ParProv_Song_2025PathMammal_BMR, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammalTree1$Ainv))
summary(MammalMCMC_BMR)


MammalMCMC_BMQ<-MCMCglmm(  fixed = Brain~Body  +Newborn + BMQ + Tb +Ta,
                           random = ~Species , data = ParProv_Song_2025PathMammal_BMR, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammalTree1$Ainv))
summary(MammalMCMC_BMQ)



##bird

ParProv_Song_2025PathBird_BMR<-ParProv_Song_2025[which(ParProv_Song_2025$Class %in% c( "Birds")
                                                       & !is.na(ParProv_Song_2025$Tb1) 
                                                       & !is.na(ParProv_Song_2025$Residual_bmr) 
                                                       & !is.na(ParProv_Song_2025$Temp) ),
                                                 c( "Species", "Residual_BS",
                                                    "Residual_OS",
                                                    "Bodymass.BdM..g.",
                                                    "Residual_bmr", 
                                                    "Tb","Temp" 
                                                 )]
nrow(ParProv_Song_2025PathBird_BMR)
summary(ParProv_Song_2025PathBird_BMR)


ParProv_Song_2025PathBird_BMR$Brain <-scale((ParProv_Song_2025PathBird_BMR$Residual_BS), center = T)
ParProv_Song_2025PathBird_BMR$Newborn<-scale((ParProv_Song_2025PathBird_BMR$Residual_OS), center = T)
ParProv_Song_2025PathBird_BMR$Body<-scale(log10(ParProv_Song_2025PathBird_BMR$Bodymass.BdM..g.), center = T)
ParProv_Song_2025PathBird_BMR$BMR<-scale((ParProv_Song_2025PathBird_BMR$Residual_bmr), center = T)
ParProv_Song_2025PathBird_BMR$Tb<-scale((ParProv_Song_2025PathBird_BMR$Tb), center = T)
ParProv_Song_2025PathBird_BMR$Ta<-scale((ParProv_Song_2025PathBird_BMR$Temp), center = T)
summary(ParProv_Song_2025PathBird_BMR)

Birdtree
BirdTree1 <- drop.tip(Birdtree , Birdtree$tip.label[which(Birdtree$tip.label %ni% rownames(ParProv_Song_2025PathBird_BMR))])
BirdTree1

inv.phyloBirdTree1 <- inverseA(chronoMPL(BirdTree1) ,nodes="TIPS",scale=TRUE)

BirdMCMC_BMR<-MCMCglmm(  fixed = Brain~Body  +Newborn + BMR + Tb +Ta,
                         random = ~Species , data = ParProv_Song_2025PathBird_BMR, prior = my.prior, 
                         nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree1$Ainv))
summary(BirdMCMC_BMR)


# Path analysis
m_bmr4 <- define_model_set(
  ##1##
  # BMR -> Tb
  m1_1  = c( Brain ~ Body + Newborn +BMR + Tb + Ta , 
             BMR ~   Body + Ta  ,
             Tb ~  Ta+ Body +BMR,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  
  # BMR ! - > Tb  
  m2_1  = c( Brain ~ Body + Newborn + BMR + Tb + Ta , 
             BMR ~   Body + Ta  ,
             Tb ~ Ta+Body ,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR - > Tb ; BMR !- > Brain 
  m3_1  = c( Brain ~ Body + Newborn + Tb + Ta , 
             BMR ~   Body + Ta  ,
             Tb ~ Ta+Body +BMR,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR !- > Tb ; BMR !- > Brain 
  m4_1  = c( Brain ~ Body + Newborn + Tb +Ta , 
             BMR ~   Body + Ta  ,
             Tb ~Ta+ Body ,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR - > Tb ; Tb !- > Brain 
  m5_1  = c( Brain ~ Body + Newborn + BMR + Ta , 
             BMR ~   Body + Ta  ,
             Tb ~ Ta+Body +BMR,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR !- > Tb ; Tb !- > Brain 
  m6_1  = c( Brain ~ Body + Newborn + Tb +Ta , 
             BMR ~   Body + Ta  ,
             Tb ~Ta+ Body ,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR - > Tb ; BMR !- > Brain ; Ta !- > Brain 
  m7_1  = c( Brain ~ Body + Newborn + Tb  , 
             BMR ~   Body + Ta  ,
             Tb ~Ta+ Body +BMR,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR !- > Tb ; BMR !- > Brain ; Ta !- > Brain 
  m8_1  = c( Brain ~ Body+ Newborn  + Tb  , 
             BMR ~   Body + Ta  ,
             Tb ~Ta+ Body ,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ) ,
  
  # BMR - > Tb ; Tb !- > Brain ; Ta !- > Brain 
  m9_1  = c( Brain ~ Body + Newborn + BMR  , 
             BMR ~   Body + Ta  ,
             Tb ~Ta+ Body +BMR,
             Newborn~ Body + Tb+ Ta ,
             Body ~ Ta
  ),
  # BMR !- > Tb ; Tb !- > Brain ; Ta !- > Brain 
  m10_1  = c( Brain ~ Body+ Newborn  + BMR  , 
              BMR ~   Body + Ta  ,
              Tb ~Ta+ Body ,
              Newborn ~ Body + Tb+ Ta ,
              Body ~ Ta
  ) 
  
  
)

#Mammal Path
nrow(ParProv_Song_2025PathMammal_BMR) # 240
p_mammal<-NA
p_mammal <- phylo_path(m_bmr4, ParProv_Song_2025PathMammal_BMR, MammalTree); p_mammal
s <- summary(p_mammal); s
#write.csv(s, paste0(PathLocation, "Mammal_Path_model.csv"))
plot(s)

b <- best(p_mammal)
plot(b)
b

avg_mammal <- average(p_mammal, cut_off = 2, avg_method = "conditional")
avg_mammal


#Bird
summary(ParProv_Song_2025PathBird_BMR)
nrow(ParProv_Song_2025PathBird_BMR)
p_bird<-NA
p_bird <- phylo_path(m_bmr4, ParProv_Song_2025PathBird_BMR, Birdtree); p_bird

s <- summary(p_bird); s
#write.csv(s, paste0(PathLocation,"Bird_Path_model.csv"))
plot(s)

b <- best(p_bird)
plot(b)
b

avg_bird <- average(p_bird, cut_off = 2, avg_method = "conditional")
avg_bird
plot(avg_bird, algorithm = 'mds', curvature = 0.1)

#####


##################################################
######Residual values results Table S11-S15########
##################################################

#prepare data
#Lampreys
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Lampreys" ),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Lampreys" ),] )
summary(a)
ParProv_Song_2025[which(ParProv_Song_2025$Class =="Lampreys"  ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which(ParProv_Song_2025$Class =="Lampreys" ),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Lampreys" ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Lampreys" ),]$Residual_OS<-a$residuals

#Cartilaginous fishes
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Cartilaginous fishes" ),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Cartilaginous fishes" ),] )
summary(a)
ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes" ),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Cartilaginous fishes" ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Cartilaginous fishes" ),]$Residual_OS<-a$residuals

#Ray-finned fishes
ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes" ),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes"  ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes"   ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes" ),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes"  ),]  )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class == "Ray-finned fishes" ),] $Residual_OS<-a$residuals

#Amphibians
ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians" ),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians" ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Amphibians"   ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which( ParProv_Song_2025$Class =="Amphibians" ),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class =="Amphibians"  ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class =="Amphibians" ),]$Residual_OS<-a$residuals

 
#Reptiles

ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles" ),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"  ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"  ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles" ),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass ),ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"   ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"  ),]$Residual_OS<-a$residuals

#Birds
ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Birds"  ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class=="Birds" ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which(ParProv_Song_2025$Class=="Birds"),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Birds"   ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class=="Birds"   ),]$Residual_OS<-a$residuals


#Mammals

ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$Residual_BS<-NA
a<-NA
a<-lm(log10(BrainMass.BrM..g.)~log10(Bodymass.BdM..g.),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Mammals" ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class=="Mammals" ),]$Residual_BS<-a$residuals

ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals"),]$Residual_OS<-NA
a<-NA
a<-lm(log10(Offspring)~log10(AdultMass),ParProv_Song_2025[which( ParProv_Song_2025$Class=="Mammals"  ),] )
summary(a)
ParProv_Song_2025[which( ParProv_Song_2025$Class=="Mammals"  ),]$Residual_OS<-a$residuals


##Table S11####

#Lampreys
# Brain size ~ Offspring size 
inv.phylojawlessTree <- inverseA(chronoMPL(jawlessTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025LampreysMCMM_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Lampreys"  )  ,]
nrow(ParProv_Song_2025LampreysMCMM_EA) # 14
ParProv_Song_2025LampreysMCMM_EA$Residual_BS<-scale(ParProv_Song_2025LampreysMCMM_EA$Residual_BS , center = T)
ParProv_Song_2025LampreysMCMM_EA$Residual_OS<-scale(ParProv_Song_2025LampreysMCMM_EA$Residual_OS , center = T)
ParProv_Song_2025LampreysMCMM_EA$Body<-NA
ParProv_Song_2025LampreysMCMM_EA$Body<-scale(log10(ParProv_Song_2025LampreysMCMM_EA$Bodymass.BdM..g.) , center = T) 

LampreyMCMCALL<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                           random = ~Species , data = ParProv_Song_2025LampreysMCMM_EA, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylojawlessTree$Ainv))
summary(LampreyMCMCALL)

#Cartilaginous fishes
# Brain size ~ Offspring size 

ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes" )  ,]
nrow(ParProv_Song_2025SharkMCMC) # 120
ParProv_Song_2025SharkMCMC$Residual_BS<-scale(ParProv_Song_2025SharkMCMC$Residual_BS , center = T)
ParProv_Song_2025SharkMCMC$Residual_OS<-scale(ParProv_Song_2025SharkMCMC$Residual_OS , center = T)
ParProv_Song_2025SharkMCMC$Body<-NA
ParProv_Song_2025SharkMCMC$Body<-scale(log10(ParProv_Song_2025SharkMCMC$Bodymass.BdM..g.) , center = T)

inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)

sharkMCMCAll<-MCMCglmm(  fixed = Residual_BS~ Body  +Residual_OS  ,
                         random = ~Species , data = ParProv_Song_2025SharkMCMC, prior = my.prior, 
                         nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCAll)

#Ray-finned fishes
# Brain size ~ Offspring size 

inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes" )  ,]
nrow(ParProv_Song_2025fishMCMC) # 210
ParProv_Song_2025fishMCMC$Residual_BS<-scale(ParProv_Song_2025fishMCMC$Residual_BS , center = T)
ParProv_Song_2025fishMCMC$Residual_OS<-scale(ParProv_Song_2025fishMCMC$Residual_OS , center = T)
ParProv_Song_2025fishMCMC$Body<-NA
ParProv_Song_2025fishMCMC$Body<-scale(log10(ParProv_Song_2025fishMCMC$Bodymass.BdM..g.) , center = T)

fishMCMCAll<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                        random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMCAll)

#Amphibians
# Brain size ~ Offspring size 
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025AmphibainMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC) #  130
ParProv_Song_2025AmphibainMCMC$Residual_BS<-scale(ParProv_Song_2025AmphibainMCMC$Residual_BS , center = T)
ParProv_Song_2025AmphibainMCMC$Residual_OS<-scale(ParProv_Song_2025AmphibainMCMC$Residual_OS , center = T)
ParProv_Song_2025AmphibainMCMC$Body<-NA
ParProv_Song_2025AmphibainMCMC$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC$Bodymass.BdM..g.) , center = T)

AmphibainMCMCAll<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                             random = ~Species , data = ParProv_Song_2025AmphibainMCMC, prior = my.prior, 
                             nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMCAll)


#Reptiles
# Brain size ~ Offspring size 
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Reptiles" & !is.na(ParProv_Song_2025$Residual_BS)   )  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 170
ParProv_Song_2025ReptileMCMC$Residual_BS<-scale(ParProv_Song_2025ReptileMCMC$Residual_BS , center = T)
ParProv_Song_2025ReptileMCMC$Residual_OS<-scale(ParProv_Song_2025ReptileMCMC$Residual_OS , center = T)
ParProv_Song_2025ReptileMCMC$Body<-NA
ParProv_Song_2025ReptileMCMC$Body<-scale(log10(ParProv_Song_2025ReptileMCMC$Bodymass.BdM..g.) , center = T)

ReptileMCMCAll<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                           random = ~Species , data = ParProv_Song_2025ReptileMCMC, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMCAll)


#Birds
# Brain size ~ Offspring size 
inv.phyloBirdTree <- inverseA(chronoMPL(Birdtree) ,nodes="TIPS",scale=TRUE)

ParProvBirdAllMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Birds"    )  ,]
nrow(ParProvBirdAllMCMC) # 1136
ParProvBirdAllMCMC$Residual_BS<-scale(ParProvBirdAllMCMC$Residual_BS , center = T)
ParProvBirdAllMCMC$Residual_OS<-scale(ParProvBirdAllMCMC$Residual_OS , center = T)
ParProvBirdAllMCMC$Body<-NA
ParProvBirdAllMCMC$Body<-scale(log10(ParProvBirdAllMCMC$Bodymass.BdM..g.) , center = T)

BirdMCMCAll<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                        random = ~Species , data = ParProvBirdAllMCMC, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree$Ainv))
summary(BirdMCMCAll)

#Mammals
# Brain size ~ Offspring size 
inv.phyloMammTree <- inverseA(chronoMPL(MammalTree) ,nodes="TIPS",scale=TRUE)

ParProv_Song_2025Mamm<-ParProv_Song_2025[which(ParProv_Song_2025$Class=="Mammals" ),]
nrow(ParProv_Song_2025Mamm) #  855
ParProv_Song_2025Mamm$Residual_BS<-scale(ParProv_Song_2025Mamm$Residual_BS , center = T)
ParProv_Song_2025Mamm$Residual_OS<-scale(ParProv_Song_2025Mamm$Residual_OS , center = T)
ParProv_Song_2025Mamm$Body<-NA
ParProv_Song_2025Mamm$Body<-scale(log10(ParProv_Song_2025Mamm$Bodymass.BdM..g.) , center = T)

MammMCMC855<-MCMCglmm(  fixed = Residual_BS~Body +Residual_OS  ,
                        random = ~Species , data = ParProv_Song_2025Mamm, prior = my.prior, 
                        nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammTree$Ainv))

summary(MammMCMC855)

#####


##Table S12####
#cacualted residual brain and new-born size across all vertebrates
a<- phylolm( log10(BrainMass.BrM..g.)  ~    log10(Bodymass.BdM..g.)    , ParProv_Song_2025  , phy = AllVerTree, model="lambda" )
ParProv_Song_2025$Residual_BS_all<-NA
for(i in 1: nrow(ParProv_Song_2025)){
  ParProv_Song_2025[i,]$Residual_BS_all<- a$residuals[which(names(a$residuals) == ParProv_Song_2025[i,]$Species)]
}

a<- phylolm( log10(Offspring)  ~    log10(AdultMass)    , ParProv_Song_2025  , phy = AllVerTree, model="lambda" )
ParProv_Song_2025$Residual_OS_all<-NA
for(i in 1: nrow(ParProv_Song_2025)){
  ParProv_Song_2025[i,]$Residual_OS_all<- a$residuals[which(names(a$residuals) == ParProv_Song_2025[i,]$Species)]
}

##all vertebrates

ParProv_Song_2025_All<-ParProv_Song_2025[which(!is.na(ParProv_Song_2025$Tb) ), 
                                          c("Species", "Class", "AdultMass", "Residual_BS_all", "Residual_OS_all", "Tb")]
nrow(ParProv_Song_2025_All) #1059
ParProv_Song_2025_All$AdultMass<-scale(log10(ParProv_Song_2025_All$AdultMass), center = T)
ParProv_Song_2025_All$Residual_BS_all<-scale(ParProv_Song_2025_All$Residual_BS_all, center = T)
ParProv_Song_2025_All$Residual_OS_all<-scale( ParProv_Song_2025_All$Residual_OS_all, center = T)
ParProv_Song_2025_All$Tb<-scale(ParProv_Song_2025_All$Tb, center = T)
summary(ParProv_Song_2025_All)

inv.phyloAllVerTree <- inverseA(chronoMPL(AllVerTree) ,nodes="TIPS",scale=TRUE)

AllVerMCMC1059<-MCMCglmm(  fixed = Residual_BS_all~  AdultMass +Residual_OS_all*Tb    ,
                           random = ~Species , data = ParProv_Song_2025_All, prior = my.prior, 
                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAllVerTree$Ainv))

summary(AllVerMCMC1059)


#####


##Table S13####
#Cartilaginous fishes
inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025SharkMCMC) # 120

sharkMCMCOff_Care<-MCMCglmm(  fixed = Residual_OS~ log10(AdultMass)  +EggCare5Cat  ,
                              random = ~Species , data = ParProv_Song_2025SharkMCMC, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCOff_Care)



sharkMCMCOff_Care_LB_Prov<-MCMCglmm(  fixed = Residual_OS~ log10(AdultMass)  +EggCare5Cat  ,
                                      random = ~Species , data = ParProv_Song_2025SharkMCMC[which(ParProv_Song_2025SharkMCMC$EggCare5Cat%in% c("3Bear", "4Pre-hat. prov.")),], prior = my.prior, 
                                      nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))

summary(sharkMCMCOff_Care_LB_Prov)

#Ray-finned fishes
inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes" & !is.na(ParProv_Song_2025$EggCare5Cat) )  ,]
nrow(ParProv_Song_2025fishMCMC) # 196

fishMCMC_Off_Care<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                              random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMC_Off_Care)



fishMCMC_Off_Care_EG_LB<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                                    random = ~Species , data = ParProv_Song_2025fishMCMC[which(ParProv_Song_2025fishMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),], prior = my.prior, 
                                    nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(fishMCMC_Off_Care_EG_LB)


#Amphibians
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)
#Offspring size ~ Egg care
ParProv_Song_2025AmphibainMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC) #  119

AmphibainMCMC_Off_care<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                                   random = ~Species , data = ParProv_Song_2025AmphibainMCMC, prior = my.prior, 
                                   nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_Off_care)


AmphibainMCMC_Off_care_EG_LB<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                                         random = ~Species , data = ParProv_Song_2025AmphibainMCMC[which(ParProv_Song_2025AmphibainMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),], prior = my.prior, 
                                         nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_Off_care_EG_LB)

#Reptiles
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)
# Offspring~Egg care
ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Reptiles"    & !is.na(ParProv_Song_2025$EggCare5Cat)  )  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 117

ReptileMCMC_Off_Care<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                                 random = ~Species , data = ParProv_Song_2025ReptileMCMC  , prior = my.prior, 
                                 nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))

summary(ReptileMCMC_Off_Care)

ReptileMCMC_Off_Care_EG_LB<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +EggCare5Cat  ,
                                       random = ~Species , data = ParProv_Song_2025ReptileMCMC[which(ParProv_Song_2025ReptileMCMC$EggCare5Cat %in% c("2Guard", "3Bear")),]  , prior = my.prior, 
                                       nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))

summary(ReptileMCMC_Off_Care_EG_LB)

#####


##Table S14####
#Cartilaginous fishes
# Brain size ~ Offspring size 
inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)

# Abandon
ParProv_Song_2025SharkMCMM_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes" & ParProv_Song_2025$EggCare5Cat =="1Abandon" )  ,]
nrow(ParProv_Song_2025SharkMCMM_EA) # 26
ParProv_Song_2025SharkMCMM_EA$Residual_BS<-scale(ParProv_Song_2025SharkMCMM_EA$Residual_BS , center = T)
ParProv_Song_2025SharkMCMM_EA$Residual_OS<-scale(ParProv_Song_2025SharkMCMM_EA$Residual_OS , center = T)
ParProv_Song_2025SharkMCMM_EA$Body<-NA
ParProv_Song_2025SharkMCMM_EA$Body<-scale(log10(ParProv_Song_2025SharkMCMM_EA$Bodymass.BdM..g.) , center = T) 

sharkMCMC_EA_scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                               random = ~Species , data = ParProv_Song_2025SharkMCMM_EA, prior = my.prior, 
                               nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_EA_scale)

#Bear
ParProv_Song_2025SharkMCMM_LB<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"   & ParProv_Song_2025$EggCare5Cat =="3Bear" )  ,]
nrow(ParProv_Song_2025SharkMCMM_LB) # 36
ParProv_Song_2025SharkMCMM_LB$Residual_BS<-scale(ParProv_Song_2025SharkMCMM_LB$Residual_BS , center = T)
ParProv_Song_2025SharkMCMM_LB$Residual_OS<-scale(ParProv_Song_2025SharkMCMM_LB$Residual_OS , center = T)
ParProv_Song_2025SharkMCMM_LB$Body<-NA
ParProv_Song_2025SharkMCMM_LB$Body<-scale(log10(ParProv_Song_2025SharkMCMM_LB$Bodymass.BdM..g.) , center = T) 

sharkMCMC_Bear_scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                 random = ~Species , data = ParProv_Song_2025SharkMCMM_LB, prior = my.prior, 
                                 nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_Bear_scale)

#Pre-hatching provision
ParProv_Song_2025SharkMCMM_PreH<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & ParProv_Song_2025$EggCare5Cat =="4Pre-hat. prov." )  ,]
nrow(ParProv_Song_2025SharkMCMM_PreH) # 58
ParProv_Song_2025SharkMCMM_PreH$Residual_BS<-scale(ParProv_Song_2025SharkMCMM_PreH$Residual_BS , center = T)
ParProv_Song_2025SharkMCMM_PreH$Residual_OS<-scale(ParProv_Song_2025SharkMCMM_PreH$Residual_OS , center = T)
ParProv_Song_2025SharkMCMM_PreH$Body<-NA
ParProv_Song_2025SharkMCMM_PreH$Body<-scale(log10(ParProv_Song_2025SharkMCMM_PreH$Bodymass.BdM..g.) , center = T) 

sharkMCMCAllPre_H_scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                    random = ~Species , data = ParProv_Song_2025SharkMCMM_PreH, prior = my.prior, 
                                    nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMCAllPre_H_scale)


#Ray-finned fishe
# Brain size ~ Offspring size 
inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025fishMCMC_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & ParProv_Song_2025$EggCare5Cat =="1Abandon"  )  ,]
nrow(ParProv_Song_2025fishMCMC_EA) # 107
ParProv_Song_2025fishMCMC_EA$Residual_BS<-scale(ParProv_Song_2025fishMCMC_EA$Residual_BS , center = T)
ParProv_Song_2025fishMCMC_EA$Residual_OS<-scale(ParProv_Song_2025fishMCMC_EA$Residual_OS , center = T)
ParProv_Song_2025fishMCMC_EA$Body<-NA
ParProv_Song_2025fishMCMC_EA$Body<-scale(log10(ParProv_Song_2025fishMCMC_EA$Bodymass.BdM..g.) , center = T)

fishMCMC_EA_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                              random = ~Species , data = ParProv_Song_2025fishMCMC_EA, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_EA_Scale)

#Guard
ParProv_Song_2025fishMCMC_EG<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & ParProv_Song_2025$EggCare5Cat =="2Guard" )  ,]
nrow(ParProv_Song_2025fishMCMC_EG) # 42
ParProv_Song_2025fishMCMC_EG$Residual_BS<-scale(ParProv_Song_2025fishMCMC_EG$Residual_BS , center = T)
ParProv_Song_2025fishMCMC_EG$Residual_OS<-scale(ParProv_Song_2025fishMCMC_EG$Residual_OS , center = T)
ParProv_Song_2025fishMCMC_EG$Body<-NA
ParProv_Song_2025fishMCMC_EG$Body<-scale(log10(ParProv_Song_2025fishMCMC_EG$Bodymass.BdM..g.) , center = T)

fishMCMC_EG_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                              random = ~Species , data = ParProv_Song_2025fishMCMC_EG, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_EG_Scale)

#Bear
ParProv_Song_2025fishMCMC_LB<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"   & ParProv_Song_2025$EggCare5Cat =="3Bear"  )  ,]
nrow(ParProv_Song_2025fishMCMC_LB) # 47
ParProv_Song_2025fishMCMC_LB$Residual_BS<-scale(ParProv_Song_2025fishMCMC_LB$Residual_BS , center = T)
ParProv_Song_2025fishMCMC_LB$Residual_OS<-scale(ParProv_Song_2025fishMCMC_LB$Residual_OS , center = T)
ParProv_Song_2025fishMCMC_LB$Body<-NA
ParProv_Song_2025fishMCMC_LB$Body<-scale(log10(ParProv_Song_2025fishMCMC_LB$Bodymass.BdM..g.) , center = T)

fishMCMC_LB_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                              random = ~Species , data = ParProv_Song_2025fishMCMC_LB, prior = my.prior, 
                              nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMC_LB_Scale)

#Amphibians
# Brain size ~ Offspring size
inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025AmphibainMCMC_EA<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"   & ParProv_Song_2025$EggCare5Cat == "1Abandon"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC_EA) #  90
ParProv_Song_2025AmphibainMCMC_EA$Residual_BS<-scale(ParProv_Song_2025AmphibainMCMC_EA$Residual_BS , center = T)
ParProv_Song_2025AmphibainMCMC_EA$Residual_OS<-scale(ParProv_Song_2025AmphibainMCMC_EA$Residual_OS , center = T)
ParProv_Song_2025AmphibainMCMC_EA$Body<-NA
ParProv_Song_2025AmphibainMCMC_EA$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC_EA$Bodymass.BdM..g.) , center = T)

AmphibainMCMC_EA_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                   random = ~Species , data = ParProv_Song_2025AmphibainMCMC_EA, prior = my.prior, 
                                   nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_EA_Scale)

#Guard
ParProv_Song_2025AmphibainMCMC_EG<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Amphibians"  & ParProv_Song_2025$EggCare5Cat == "2Guard"  )  ,]
nrow(ParProv_Song_2025AmphibainMCMC_EG) #  25
ParProv_Song_2025AmphibainMCMC_EG$Residual_BS<-scale(ParProv_Song_2025AmphibainMCMC_EG$Residual_BS , center = T)
ParProv_Song_2025AmphibainMCMC_EG$Residual_OS<-scale(ParProv_Song_2025AmphibainMCMC_EG$Residual_OS , center = T)
ParProv_Song_2025AmphibainMCMC_EG$Body<-NA
ParProv_Song_2025AmphibainMCMC_EG$Body<-scale(log10(ParProv_Song_2025AmphibainMCMC_EG$Bodymass.BdM..g.) , center = T)

AmphibainMCMC_EG_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                   random = ~Species , data = ParProv_Song_2025AmphibainMCMC_EG, prior = my.prior, 
                                   nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))

summary(AmphibainMCMC_EG_Scale)


#Reptiles
#Brain size ~ Offspring size
inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)

#Abandon
ParProv_Song_2025ReptileMCMC_EA<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"    & ParProv_Song_2025$EggCare5Cat=="1Abandon" )  ,]
nrow(ParProv_Song_2025ReptileMCMC_EA) # 66
ParProv_Song_2025ReptileMCMC_EA$Residual_BS<-scale(ParProv_Song_2025ReptileMCMC_EA$Residual_BS , center = T)
ParProv_Song_2025ReptileMCMC_EA$Residual_OS<-scale(ParProv_Song_2025ReptileMCMC_EA$Residual_OS , center = T)
ParProv_Song_2025ReptileMCMC_EA$Body<-NA
ParProv_Song_2025ReptileMCMC_EA$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_EA$Bodymass.BdM..g.) , center = T)


ReptileMCMC_EA_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                 random = ~Species , data = ParProv_Song_2025ReptileMCMC_EA, prior = my.prior, 
                                 nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_EA_Scale)


#2Guard
ParProv_Song_2025ReptileMCMC_EG<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"    & ParProv_Song_2025$EggCare5Cat =="2Guard"   )  ,]
nrow(ParProv_Song_2025ReptileMCMC_EG) # 21
ParProv_Song_2025ReptileMCMC_EG$Residual_BS<-scale(ParProv_Song_2025ReptileMCMC_EG$Residual_BS , center = T)
ParProv_Song_2025ReptileMCMC_EG$Residual_OS<-scale(ParProv_Song_2025ReptileMCMC_EG$Residual_OS , center = T)
ParProv_Song_2025ReptileMCMC_EG$Body<-NA
ParProv_Song_2025ReptileMCMC_EG$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_EG$Bodymass.BdM..g.) , center = T)

ReptileMCMC_EG_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                 random = ~Species , data = ParProv_Song_2025ReptileMCMC_EG, prior = my.prior, 
                                 nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_EG_Scale)

#Bear
ParProv_Song_2025ReptileMCMC_LB<-ParProv_Song_2025[which( ParProv_Song_2025$Class   =="Reptiles"   &  ParProv_Song_2025$EggCare5Cat=="3Bear"  )  ,]
nrow(ParProv_Song_2025ReptileMCMC_LB) # 34
ParProv_Song_2025ReptileMCMC_LB$Residual_BS<-scale(ParProv_Song_2025ReptileMCMC_LB$Residual_BS , center = T)
ParProv_Song_2025ReptileMCMC_LB$Residual_OS<-scale(ParProv_Song_2025ReptileMCMC_LB$Residual_OS , center = T)
ParProv_Song_2025ReptileMCMC_LB$Body<-NA
ParProv_Song_2025ReptileMCMC_LB$Body<-scale(log10(ParProv_Song_2025ReptileMCMC_LB$Bodymass.BdM..g.) , center = T)

ReptileMCMC_LB_Scale<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                 random = ~Species , data = ParProv_Song_2025ReptileMCMC_LB, prior = my.prior, 
                                 nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMC_LB_Scale)

#Birds
# Brain size ~ Offspring size
# all incubation
inv.phyloBirdTree <- inverseA(chronoMPL(Birdtree) ,nodes="TIPS",scale=TRUE)

ParProvBirdAllMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Birds"    )  ,]
nrow(ParProvBirdAllMCMC) # 1136
ParProvBirdAllMCMC$Residual_BS<-scale(ParProvBirdAllMCMC$Residual_BS , center = T)
ParProvBirdAllMCMC$Residual_OS<-scale(ParProvBirdAllMCMC$Residual_OS , center = T)
ParProvBirdAllMCMC$Body<-NA
ParProvBirdAllMCMC$Body<-scale(log10(ParProvBirdAllMCMC$Bodymass.BdM..g.) , center = T)

BirdMCMC_Pre_Hatch_incubte<-MCMCglmm(  fixed = Residual_BS~Body  +Residual_OS  ,
                                       random = ~Species , data = ParProvBirdAllMCMC[which( ParProvBirdAllMCMC$EggCare5Cat!="1Abandon"),], prior = my.prior, 
                                       nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree$Ainv))
summary(BirdMCMC_Pre_Hatch_incubte)

#####


##Table S15####

# shark
ParProv_Song_2025SharkMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class =="Cartilaginous fishes"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025SharkMCMC) # 118
ParProv_Song_2025SharkMCMC$Residual_BS<-scale(ParProv_Song_2025SharkMCMC$Residual_BS , center = T)
ParProv_Song_2025SharkMCMC$Residual_OS<-scale(ParProv_Song_2025SharkMCMC$Residual_OS , center = T)
ParProv_Song_2025SharkMCMC$Body<-NA
ParProv_Song_2025SharkMCMC$Body<-scale(log10(ParProv_Song_2025SharkMCMC$Bodymass.BdM..g.) , center = T)

ParProv_Song_2025SharkMCMC$AdultMass<-scale(log10(ParProv_Song_2025SharkMCMC$AdultMass) , center = T)
ParProv_Song_2025SharkMCMC$Tb<-scale(ParProv_Song_2025SharkMCMC$Tb , center = T)


inv.phyloSharkTree <- inverseA(chronoMPL(SharkTree) ,nodes="TIPS",scale=TRUE)


sharkMCMC_BS_Tb<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS+ Tb  ,
                            random = ~Species , data = ParProv_Song_2025SharkMCMC[which(!is.na(ParProv_Song_2025SharkMCMC$Tb)),], prior = my.prior, 
                            nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloSharkTree$Ainv))
summary(sharkMCMC_BS_Tb)


## fish 
ParProv_Song_2025fishMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Ray-finned fishes"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025fishMCMC) # 159
ParProv_Song_2025fishMCMC$Residual_BS<-scale(ParProv_Song_2025fishMCMC$Residual_BS , center = T)
ParProv_Song_2025fishMCMC$Residual_OS<-scale(ParProv_Song_2025fishMCMC$Residual_OS , center = T)

ParProv_Song_2025fishMCMC$Body<-NA
ParProv_Song_2025fishMCMC$Body<-scale(log10(ParProv_Song_2025fishMCMC$Bodymass.BdM..g.) , center = T)

ParProv_Song_2025fishMCMC$AdultMass<-scale(log10(ParProv_Song_2025fishMCMC$AdultMass) , center = T)
ParProv_Song_2025fishMCMC$Tb<-scale(ParProv_Song_2025fishMCMC$Tb , center = T)

inv.phylofishTree <- inverseA(chronoMPL(fishtree) ,nodes="TIPS",scale=TRUE)

fishMCMCBS_Tb<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS+ Tb  , 
                          random = ~Species , data = ParProv_Song_2025fishMCMC, prior = my.prior, 
                          nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))
summary(fishMCMCBS_Tb)




## Amphibias 
ParProv_Song_2025AmphibiansMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Amphibians"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025AmphibiansMCMC) # 41
ParProv_Song_2025AmphibiansMCMC$Residual_BS<-scale(ParProv_Song_2025AmphibiansMCMC$Residual_BS , center = T)
ParProv_Song_2025AmphibiansMCMC$Residual_OS<-scale(ParProv_Song_2025AmphibiansMCMC$Residual_OS , center = T)

ParProv_Song_2025AmphibiansMCMC$Body<-NA
ParProv_Song_2025AmphibiansMCMC$Body<-scale(log10(ParProv_Song_2025AmphibiansMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025AmphibiansMCMC$AdultMass<-scale(log10(ParProv_Song_2025AmphibiansMCMC$AdultMass) , center = T)

ParProv_Song_2025AmphibiansMCMC$Tb<-scale(ParProv_Song_2025AmphibiansMCMC$Tb , center = T)

inv.phyloAmphibainTree <- inverseA(chronoMPL(AmphibiaTree) ,nodes="TIPS",scale=TRUE)

AmphibianMCMCBrain_Temp<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS +Tb  ,
                                                random = ~Species , data = ParProv_Song_2025AmphibiansMCMC, prior = my.prior, 
                                                nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloAmphibainTree$Ainv))
summary(AmphibianMCMCBrain_Temp)




## Reptile 
ParProv_Song_2025ReptileMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Reptiles"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025ReptileMCMC) # 112
ParProv_Song_2025ReptileMCMC$Residual_BS<-scale(ParProv_Song_2025ReptileMCMC$Residual_BS , center = T)
ParProv_Song_2025ReptileMCMC$Residual_OS<-scale(ParProv_Song_2025ReptileMCMC$Residual_OS , center = T)

ParProv_Song_2025ReptileMCMC$Body<-NA
ParProv_Song_2025ReptileMCMC$Body<-scale(log10(ParProv_Song_2025ReptileMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025ReptileMCMC$AdultMass<-scale(log10(ParProv_Song_2025ReptileMCMC$AdultMass) , center = T)

ParProv_Song_2025ReptileMCMC$Tb<-scale(ParProv_Song_2025ReptileMCMC$Tb , center = T)

inv.phyloReptileTree <- inverseA(chronoMPL(CrocoTurtleMCCTree) ,nodes="TIPS",scale=TRUE)


ReptileMCMCBrain_Temp<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS +Tb  ,
                                    random = ~Species , data = ParProv_Song_2025ReptileMCMC, prior = my.prior, 
                                    nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloReptileTree$Ainv))
summary(ReptileMCMCBrain_Temp)

## Birds 
ParProv_Song_2025BirdMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Birds"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025BirdMCMC) # 348
ParProv_Song_2025BirdMCMC$Residual_BS<-scale(ParProv_Song_2025BirdMCMC$Residual_BS , center = T)
ParProv_Song_2025BirdMCMC$Residual_OS<-scale(ParProv_Song_2025BirdMCMC$Residual_OS , center = T)

ParProv_Song_2025BirdMCMC$Body<-NA
ParProv_Song_2025BirdMCMC$Body<-scale(log10(ParProv_Song_2025BirdMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025BirdMCMC$AdultMass<-scale(log10(ParProv_Song_2025BirdMCMC$AdultMass) , center = T)

ParProv_Song_2025BirdMCMC$Tb<-scale(ParProv_Song_2025BirdMCMC$Tb , center = T)

BirdTree4<-drop.tip(Birdtree, Birdtree$tip.label[! Birdtree$tip.label %in% ParProv_Song_2025BirdMCMC$Species],
                    root.edge = F, rooted = is.rooted(Birdtree))
inv.phyloBirdTree4 <- inverseA(chronoMPL(BirdTree4) ,nodes="TIPS",scale=TRUE)

 
BirdMCMCBrain_Temp<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS +Tb  ,
                                  random = ~Species , data = ParProv_Song_2025BirdMCMC, prior = my.prior, 
                                  nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloBirdTree4$Ainv))
summary(BirdMCMCBrain_Temp)

## Mammals 
ParProv_Song_2025MammalMCMC<-ParProv_Song_2025[which(ParProv_Song_2025$Class == "Mammals"  & !is.na(ParProv_Song_2025$Tb))  ,]
nrow(ParProv_Song_2025MammalMCMC) # 274
ParProv_Song_2025MammalMCMC$Residual_BS<-scale(ParProv_Song_2025MammalMCMC$Residual_BS , center = T)
ParProv_Song_2025MammalMCMC$Residual_OS<-scale(ParProv_Song_2025MammalMCMC$Residual_OS , center = T)

ParProv_Song_2025MammalMCMC$Body<-NA
ParProv_Song_2025MammalMCMC$Body<-scale(log10(ParProv_Song_2025MammalMCMC$Bodymass.BdM..g.) , center = T)
ParProv_Song_2025MammalMCMC$AdultMass<-scale(log10(ParProv_Song_2025MammalMCMC$AdultMass) , center = T)

ParProv_Song_2025MammalMCMC$Tb<-scale(ParProv_Song_2025MammalMCMC$Tb , center = T)


MammalTree1<-drop.tip(MammalTree, MammalTree$tip.label[! MammalTree$tip.label %in% ParProv_Song_2025MammalMCMC$Species],
                      root.edge = F, rooted = is.rooted(MammalTree))
inv.phyloMammTree1 <- inverseA(chronoMPL(MammalTree1) ,nodes="TIPS",scale=TRUE)



MammalMCMCBrain_Temp<-MCMCglmm(  fixed = Residual_BS ~  Body +Residual_OS +Tb  ,
                               random = ~Species , data = ParProv_Song_2025MammalMCMC, prior = my.prior, 
                               nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phyloMammTree1$Ainv))
summary(MammalMCMCBrain_Temp)



#####


##Table S16, Fertilization####

ParProvFishMCMCFertilization<-ParProv_Song_2025[which(ParProv_Song_2025$Class   =="Ray-finned fishes"    & !is.na(ParProv_Song_2025$Fertilization)  )  ,]
nrow(ParProvFishMCMCFertilization) # 186
table(ParProvFishMCMCFertilization$Fertilization)

#Offspring
FishMCMCFertilizationMCMC_Off<-MCMCglmm(  fixed = Residual_OS ~ log10(AdultMass)  +Fertilization  ,
                                          random = ~Species , data = ParProvFishMCMCFertilization  , prior = my.prior, 
                                          nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(FishMCMCFertilizationMCMC_Off)

#Brain
FishMCMCFertilizationMCMC_Brain<-MCMCglmm(  fixed = Residual_BS ~ log10(Bodymass.BdM..g.)  +Fertilization  ,
                                            random = ~Species , data = ParProvFishMCMCFertilization  , prior = my.prior, 
                                            nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(FishMCMCFertilizationMCMC_Brain)

#Adult mass
FishMCMCFertilizationMCMC_Body<-MCMCglmm(  fixed = log10(AdultMass) ~ Fertilization  ,
                                           random = ~Species , data = ParProvFishMCMCFertilization  , prior = my.prior, 
                                           nitt = 500000, thin = 40, burnin = 100000, verbose = F, ginverse = list(Species = inv.phylofishTree$Ainv))

summary(FishMCMCFertilizationMCMC_Body)

#####


