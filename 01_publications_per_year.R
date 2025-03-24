
# run this first
if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(dplyr, RefManageR, ggplot2, patchwork, ggpubr, rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path
# set working directory to active directory
setwd(dirname(current_path ))

getwd()

# create subfolder for saving figures
dir.create(file.path(getwd(), "figures")) 


##############################################################
# Publications per year
##############################################################

# reading from Latex .bib file of the list of EMF studies included in the review
bib <- ReadBib("JabRef_5a_insects_EMF.bib", check = "warn")
df <- as.data.frame(bib, row.names = TRUE)
df$Identifier <- row.names(df)
df$year <- substr(df$date, 1, 4)
df$year <- as.numeric(df$year)
#View(df)


p1 <- ggplot(data = df, aes(year)) + 
  labs(title = "Insects and EMF", y = "Number of publications", x = "Year") + 
  theme_classic() + geom_histogram(binwidth = 0.5, fill="red") +
  scale_x_continuous(limits = c(1980, 2023)) +
  theme(legend.position = "bottom")
p1

# reading from Latex .bib file of the list of magnetic sense studies included in the review
bib2 <- ReadBib("JabRef_5b_magnetic_sense.bib", check = "error")
df2 <- as.data.frame(bib2, row.names = TRUE)
df2$Identifier <- row.names(df2)
df2$year <- substr(df2$date, 1, 4)
df2$year <- as.numeric(df2$year)
#View(df2)

p2 <- ggplot(data = df2, aes(year)) + 
  labs(title = "Insects and magnetic sense", y = "Number of publications", x = "Year") + 
  theme_classic() + geom_histogram(binwidth = 0.5, fill="blue") +
  scale_x_continuous(limits = c(1980, 2023)) + scale_y_continuous(limits = c(0, 16)) +
  theme(legend.position = "bottom")
p2


#patchwork::align_plots(p1,p2)
patchwork = p1 + p2
patchwork[[2]] = patchwork[[2]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank() )
patchwork

#ggsave("figures/publications_all.jpg", width = 9, height = 3)

ggarrange(p1, p2, ncol = 1, nrow = 2, align = "hv")

#ggsave("figures/publications_all(vertical).jpg", width = 2.5, height = 7)


#########################################################
# Create an overlay histogram with ggplot2
#########################################################

df$type <- "EMF"
df2$type <- "Magnetic sense"
df3 <- full_join(df,df2)

ggplot(df3, aes(year, fill = type)) +
  geom_histogram(binwidth = 0.5) +
  labs(title = "", y = "Number of publications", x = "Year", fill = "Topic") +
  theme_classic() + scale_x_continuous(limits = c(1980, 2022)) +
  theme(legend.position = "bottom")

ggsave("figures/publications_overlaid.jpg", width = 3, height = 3)


