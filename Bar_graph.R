library(ggplot2)

hg19 = 28299634
hg38 = hg19

row = data.frame("sample" = "107", 
                 "wave" = "wave3",
                 "CpGs_850" = 860000,
                 "CpGs_wave" = 186000,
                 "PASS1" = 3299634,
                 "PASS10" = 299634,
                 "PASS30" = 29963,
                 "PASS50" = 2996
                )

row %>%
  ggplot( aes(y=PASS30, group=Approach, color=Approach)) +
  geom_bar() +
  ggtitle("Coverage across samples") +
  ylab("Coverage (X)") +
  xlab("Sample")


