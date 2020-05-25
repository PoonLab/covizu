x <- c(as.Date('2020-04-03'), as.Date('2020-04-07'), 
       as.Date('2020-04-10'), as.Date('2020-04-14'),
       as.Date('2020-04-19'))
y <- c(3815, 4645, 5818, 8189, 10650)
par(mar=c(5,5,1,1), bty='n')

plot(x, y, type='b', xlab='Release date', ylab='Number of genomes')


#require(xlsx)
#ack <- read.xlsx(file='~/Downloads/gisaid_cov2020_acknowledgement_table.xls',
#                 sheetIndex = 1, startRow = 3)
#ack$coldate <- as.Date(ack$Collection.date)
#write.csv(data.frame(accession=ack$Accession.ID, desc=ack$Virus.name, coldate=ack$coldate),
#          file='~/git/covizu/data/acknow.csv')

ack <- read.csv('~/git/covizu/data/acknow.csv', row.names=1)
ack$coldate <- as.Date(ack$coldate)

require(ggfree)

par(xpd=FALSE, mar=c(5,5,1,1))
plot(NA,
     xlim=c(as.Date('2020-01-01'), 
            max(ack$coldate, na.rm=T)),
     ylim=c(1, nrow(ack)),
     xlab='', ylab='Number of genomes',
     xaxt='n', bty='n', col='cadetblue', lwd=2, las=1, cex.axis=0.75)
add.grid(mode='x', bg.col = 'white', fg.col='grey90')
lines(sort(ack$coldate), 1:sum(!is.na(ack$coldate)), 
      type='s', col='cadetblue', lwd=2)


xt <- seq(as.Date('2020-01-01'), max(ack$coldate, na.rm=T), length.out=10)

axis(side=1, at=xt, labels=strftime(xt, '%b %d'), 
     las=2, cex.axis=0.7, mgp=c(3,0.6,0))
title(xlab='Collection date', line=3.)

