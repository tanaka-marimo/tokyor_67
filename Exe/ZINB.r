#--------------------------------------------------------------------------------------
#  TokyoR #67 豪雪　LT 
#
#  !! インプットファイル !!
#  ../Data/snow.csv                            : JMAHPより取得
#  ../Data/wind                                : JMAHPより取得
#
#                                                         2018.1.13 新規
#
#--------------------------------------------------------------------------------------

rm(list=ls())

library("dplyr")
library("ggplot2")
library("data.table")
library("stringr")
library("rstan")
require(shinystan)

# read
snow <- fread("../Data/snow.csv")
wind <- fread("../Data/wind.csv")

# join snow wind
df <- snow %>%
  dplyr::left_join(wind,by=c("date")) %>%
  dplyr::select(date,derection,speed,temp,nou,toukamachi,niitsu)

ymd <- str_split(df$date, pattern = "/") %>% 
  as.data.frame %>% t 
dimnames(ymd) <- NULL
ymd <- as.data.frame(ymd)
colnames(ymd) <- c("yyyy","mm","dd")
df_ymd <- cbind(df,ymd)

# NA cut, filter wnw ~ nne
fuyugata <- c("西北西", "北西", "北北西","北","北北東")
df_fuyugata <- df_ymd %>%
  dplyr::filter(derection %in% fuyugata) %>%
  na.omit

# メカニズムの想定 ヒストグラム確認;
df_fuyugata$nou %>% hist
df_fuyugata$toukamachi %>% hist
df_fuyugata$niitsu %>% hist

# SST parameter
df_fuyu_sst <- df_fuyugata %>%
  dplyr::mutate(dec = if_else(condition = mm == 12, true = 1, false = 0)) %>%
  dplyr::mutate(jan = if_else(condition = mm == 1, true = 1, false = 0)) %>%
  dplyr::mutate(feb = if_else(condition = mm == 2, true = 1, false = 0)) %>%
  dplyr::select(-date,-derection)

# Sado parameter
df_fuyu_nou <- df_fuyu_sst %>%
  dplyr::select(-niitsu,-toukamachi) %>%
  dplyr::mutate(sado=1)
names(df_fuyu_nou)[3] <- "snow"

df_fuyu_niitsu <- df_fuyu_sst %>%
  dplyr::select(-toukamachi,-nou) %>%
  dplyr::mutate(sado=0)
names(df_fuyu_niitsu)[3] <- "snow"

df_sst_sado <- rbind(df_fuyu_nou,df_fuyu_niitsu)

# stan に渡すデータ no sst ver
datastan <- list(N=nrow(df_sst_sado), Y=df_sst_sado$snow,
                 X1=df_sst_sado$temp, X2=df_sst_sado$speed,
                 S=df_sst_sado$sado)

# negative binominal ver
stanmodel <- stan_model(file='../Model/ZINB.stan')
fit <- sampling(stanmodel, data=datastan,
                pars=c('c', 'd0', 'd1', 'c_new', 'a0', 'a1', 'a2','phi','mu_new'),
                seed=1234)

# 収束診断
launch_shinystan(fit)

# 作図
imgdir <-"../Data/OUT/"
# ロジステック関数
# 準備（雪ありを１に）
pre_gg <- df_sst_sado %>%
  dplyr::mutate(umu = if_else(condition = snow > 0, true = 1, false = 0))

d0<-rstan::extract(fit)$d0%>%as.numeric() %>% mean
d1<-rstan::extract(fit)$d1%>%as.numeric() %>% mean
logist <- function(x){
  y = exp(d0 + ( d1*x)) / (1 + exp(d0 + ( d1*x)))
}

g <-ggplot(pre_gg, aes(temp, umu)) + geom_point() +
  theme_gray (base_family = "HiraKakuPro-W3") +
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) +
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) 

g <- g + stat_function(fun = logist) + xlab("羽茂の気温") +ylab("雨雪判別")
plot(g)

# nou ver
ggsave(filename=paste0(imgdir,"logst.png"),
       width=5,height=4,dpi=200)

# c_newのヒストグラム
c_new<-rstan::extract(fit)$c_new%>%as.numeric()
res<-data.frame(sample=c(c_new),pars=rep(c("佐渡島の寄与"),each=length(c_new)))
Color=c("tomato")
g <- ggplot(res,aes(x=sample,y=..density..,colour=pars,fill=pars))+
  theme_gray (base_family = "HiraKakuPro-W3") +
  geom_histogram(position="identity",alpha=0.7)+
  geom_density(position="identity",alpha=0.4)+
  scale_color_manual(values=Color)+
  scale_fill_manual(values=Color) +
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) +
  xlab("佐渡の寄与") +ylab("denstiy")

plot(g)

ggsave(filename=paste0(imgdir,"sado_nb.png"),
       width=6,height=4,dpi=200)

# 推定確認作図
mu_new <- rstan::extract(fit)$mu_new %>%as.data.frame()
mu_gg <- mu_new %>%
  dplyr::summarise_each(funs(mean),everything()) %>%
  as.numeric

res<-data.frame(sample=c(mu_gg),pars=rep(c("推定値"),each=length(mu_gg)))
Color=c("tomato")
g <- ggplot(res,aes(x=sample,y=..density..,colour=pars,fill=pars))+
  theme_gray (base_family = "HiraKakuPro-W3") +
  geom_histogram(position="identity",alpha=0.7)+
  scale_color_manual(values=Color)+
  scale_fill_manual(values=Color) +
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) +
  xlab("推定値") +ylab("denstiy")

plot(g)

ggsave(filename=paste0(imgdir,"fcst_hist.png"),
       width=5,height=4,dpi=200)

# 比較用実況値作図
res<-data.frame(sample=c(df_sst_sado$snow),pars=rep(c("実況値"),each=nrow(df_sst_sado)))
Color=c("dimgray")
g <- ggplot(res,aes(x=sample,y=..density..,colour=pars,fill=pars))+
  theme_gray (base_family = "HiraKakuPro-W3") +
  geom_histogram(position="identity",alpha=0.7)+
  scale_color_manual(values=Color)+
  scale_fill_manual(values=Color) +
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) +
  xlab("実況値") +ylab("denstiy")

plot(g)

ggsave(filename=paste0(imgdir,"obs_hist.png"),
       width=5,height=4,dpi=200)

# 比較用実況値 & 推定値　を重ねたヒストグラム作図
res<-data.frame(sample=c(df_sst_sado$snow,mu_gg),pars=rep(c("実況値","推定値"),each=nrow(df_sst_sado)))
Color=c("dimgray","tomato")
g <- ggplot(res,aes(x=sample,y=..density..,colour=pars,fill=pars))+
  theme_gray (base_family = "HiraKakuPro-W3") +
  geom_histogram(position="identity",alpha=0.7)+
  scale_color_manual(values=Color)+
  scale_fill_manual(values=Color) +
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15)) +
  xlab("降雪量") +ylab("denstiy")

plot(g)

ggsave(filename=paste0(imgdir,"obs_fcst_hist.png"),
       width=5,height=4,dpi=200)
