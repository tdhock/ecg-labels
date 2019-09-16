source("packages.R")

data(ECG, package="gfpop")

m <- function(x){
  factor(x, c("previous", "changepoint"), c(
    "Previous model", "Proposed models"))
}
arg.list <- list(
  list("b", "Q",      "down", 8000000, gap=0),
  list("Q",       "R",      "up",   0,        gap=2000),
  list("R",       "S",      "down", 0,        gap=5000),
  list("S",       "x1",     "up",   0,        gap=2000),
  list("x1",      "x2",     "up",   0,        gap=1000),
  list("x2",      "x3",     "up",   0,        gap=0),
  list("x3",      "x4",     "up",   0,        gap=0,  model=1),
  list("x4",      "y1",   "up",   0,        gap=0,  model=1),
  list("y1",    "y2",  "down", 0,        gap=0),
  list("y2",   "y3", "down", 0,        gap=0),
  list("y3",  "y4", "down", 0,        gap=0, model=1),
  list("y4",  "y5",   "down",   0,        gap=0, model=1),
  list("y5",    "b","up",   0,        gap=0))
model.dt.list <- list(
  previous=ECG$PanTompkins[, data.table(
    n.states=NA,
    states=NA,
    model=m("previous"),
    time,
    millivolts,
    state=letter
  )])
mean.dt.list <- list()
select.vec.list <- list(
  less=sapply(arg.list, function(L)is.null(L$model)),
  more=rep(TRUE, length(arg.list)))

for(states in names(select.vec.list)){
  select.vec <- select.vec.list[[states]]
  edge.list <- lapply(arg.list[select.vec], function(L){
    do.call(gfpop::Edge, L[names(L)!="model"])
  })
  prev <- edge.list[[1]]$state1
  for(edge.i in seq_along(edge.list)){
    if(edge.list[[edge.i]]$state1 != prev){
      edge.list[[edge.i]]$state1 <- prev
    }
    prev <- edge.list[[edge.i]]$state2
  }
  myGraph <- do.call(gfpop::graph, edge.list)
  fit <- gfpop::gfpop(ECG$data$millivolts, mygraph = myGraph, type = "mean")
  end.i <- fit$changepoints
  start.i <- c(1, end.i[-length(end.i)]+1)
  n.states <- length(unique(fit$states))
  segments.dt <- with(ECG$data, data.table(
    timeStart=time[start.i],
    timeEnd=time[end.i],
    state=fit$states,
    mean=fit$parameters))
  ## to show vertical lines connecting the horizontal segment means:
  mean.dt.list[[states]] <- segments.dt[, data.table(
    n.states,
    states,
    time=as.numeric(rbind(timeStart-0.5, timeEnd+0.5)),
    mean=as.numeric(rbind(mean, mean)))]
  model.dt.list[[states]] <- segments.dt[, data.table(
    n.states,
    states,
    model=m("changepoint"),
    time=ifelse(state=="Q", timeEnd, (timeStart+timeEnd)/2),
    millivolts=mean,
    state)]
}

model.dt <- do.call(rbind, model.dt.list)[state %in% c("Q", "R", "S")]
mean.dt <- do.call(rbind, mean.dt.list)
mean.dt[, States := factor(n.states, c("9", "13"))]
samples.per.second <- 250
truth.dt <- segments.dt[state=="R", list(time=(timeStart+timeEnd)/2)]
gg <- ggplot()+
  geom_vline(aes(
    xintercept=time/samples.per.second),
    color="red",
    data=truth.dt)+
  geom_text(aes(
    x, y/1e3, hjust=hjust, label="True R"),
    color="red",
    size=3,
    data=data.table(
      x=208.5, y=6500, hjust=1, label="True R", model=m("changepoint")))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(model ~ .)+
  geom_line(aes(
    time/samples.per.second, millivolts/1e3),
    color="grey50",
    data=ECG$data)+
  geom_line(aes(
    time/samples.per.second, mean/1e3, color=States, size=States),
    data=data.table(model=m("changepoint"), mean.dt),
    alpha=0.5
    )+
  scale_color_manual(values=c(
    "13"="blue",
    "9"="deepskyblue"))+
  scale_size_manual(values=c(
    "13"=0.8,
    "9"=1.5))+
  geom_label(aes(
    time/samples.per.second, millivolts/1e3,
    label=state),
    color="blue",
    size=3,
    label.padding=grid::unit(0.1, "lines"),
    alpha=0.6,
    data=model.dt)+
  coord_cartesian(xlim=c(52000, 52900)/samples.per.second, expand=FALSE)+
  xlab("Time (seconds)")+
  ylab("Electrocardiogram activity (Volts)")
png("figure-two-ecg-graphs-data.png", 8, 2.6, res=400, units="in")
print(gg)
dev.off()

## bottom part of figure
cat(tikz.tex <- data.table(myGraph)[type != "null", paste(c("
\\definecolor{deepskyblue}{RGB}{0,191,255}
\\begin{tikzpicture}[->,>=latex,shorten >=1pt,auto,node distance=1.3cm,
      thick,main node/.style={circle,draw}]
", sprintf(
  "\\node[main node, fill=%s, text=blue] (%s) %s {%s};\n",
  ifelse(select.vec.list$less, "white!50!deepskyblue", "white"),
  state1,
  ifelse(state2=="Q", "", paste0("[right of=", c(NA, state1[-.N]), "]")),
  state1), "
\\path[every node/.style={font=\\sffamily\\small}]
", sprintf(
  "(%s) edge [bend left%s] node [above] {%s} node [below] {\\scriptsize %s} (%s)\n",
  state1, ifelse(
    state2=="b", ", looseness=0.4", ""),
  sprintf(
    "$%s, %s$",
    ifelse(
      penalty==0, 0, "\\lambda"),
    ifelse(
      type=="up", "\\uparrow", "\\downarrow")),
  ifelse(parameter==0, "", parameter/1e3),
  state2), ";
\\end{tikzpicture}
"), sep="\n")], file="figure-two-ecg-graphs.tex")

