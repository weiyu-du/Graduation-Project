#pre：单细胞类型识别结果
#y：真实的单细胞类型
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
  print('F1：')
  print(f1)
  print('Mean（F1）；')
  print(mean(f1))
}
