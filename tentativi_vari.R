load("/home/simon/Documenti/GitHub/gruppoFucksia/data/dataset_fix.RData")
dim(colData(se_fix))
dim(assay(se_fix))
dim(rowData(se_fix))

install.packages(plotly)
library("plotly")

x <- 
y <- rnorm(1000)
s <- subplot(
  plot_ly(x = x, color = I("black"), type = 'histogram'), 
  plotly_empty(), 
  plot_ly(x = x, y = y, type = 'histogram2dcontour', showscale = F), 
  plot_ly(y = y, color = I("black"), type = 'histogram'),
  nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), 
  shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
)

fig <- layout(s, showlegend = FALSE)

fig

cor_fix <- abs(cor(t(assay(se_fix))))
fig <- plot_ly(z = cor_fix, type = "heatmap")

fig

tmp <- data_fix[,c(2,10:508)]
cor_tmp <- cor(tmp[,2:500])
tmp_adenoma <- data_fix[data_fix$study_condition == "adenoma",10:508]
tmp_control <- data_fix[data_fix$study_condition == "control",9:508]

cor_adenoma <- cor(tmp_adenoma)
fig <- plot_ly(z = cor_adenoma, type = "heatmap")

fig

