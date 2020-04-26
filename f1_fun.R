#pre：预测的分类结果
#y：真实的分类结果
f1_fun = function(pre,y){
  class = sort(unique(y))
  tp=NA
  fp=NA
  fn=NA
  for(i in 1:length(class)){
    tp[i] = sum(pre==class[i] & y==class[i])
    fp[i] = sum(pre==class[i] & y!=class[i])
    fn[i] = sum(pre!=class[i] & y==class[i])
  }
  f1 = 2*tp/(2*tp+fp+fn)
  names(f1) = class
  print(table(pre,y))
  print('-----------------f1----------------------')
  print(f1)
  print('--------------mean(f1)-------------------')
  print(mean(f1))
}
