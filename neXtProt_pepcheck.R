#-----------------check unique peptide
neXtProt_pepcheck <- function(Peptide,no_variant_match=T,mode='J'){
  options(timeout=9999999)
  library(httr)
  library(jsonlite)
  # 定义API的基础URL
  base_url <- "https://api.nextprot.org/"
  # 定义查询参数
  # no_variant_match <- "false" # 或 "false"，根据需求
  # mode <- "J" # 根据API文档设定
  # 构造完整的请求URL
  if(length(Peptide)>500){
    cat('Peptide is more than 500')
    sub_Peptide <- split(Peptide, rep(1:(length(Peptide) %/% 500 + 1), each = 500, length.out = length(Peptide)))
    cat(paste0('Need run ',length(sub_Peptide)))
    res <- list()
    for (i in 1:length(sub_Peptide)) {
      # if(i==355){
      #   cat(paste0('让访问暂停30秒'))
      #   Sys.sleep(30)
      # }
      request_url <- paste0(base_url, "entries/search/peptide.json?peptide=", paste0(sub_Peptide[[i]],collapse = ','), 
                            "&no-variant-match=", no_variant_match, 
                            "&mode=", mode)
      
      # 发送GET请求 添加遇到报错，重试的步骤
      max_retries <- 3
      retry_count <- 0
      while (retry_count < max_retries) {
        tryCatch({
          response <- GET(request_url)
          # 如果请求成功，跳出循环
          break
        }, error = function(e) {
          message(paste("请求失败，尝试第", retry_count + 1, "次重试"))
          retry_count <<- retry_count + 1
          # 可以在这里加入延迟，避免立即重试导致服务器压力过大
          Sys.sleep(30)  # 等待5秒后重试
        })
      }
      
      # 发送GET请求
      # response <- GET(request_url)
      
      # 检查响应状态码
      if (http_status(response)$category == "Success") {
        print(paste0('Run',i,' GET is ok'))
        # 解析JSON响应内容（假设API返回的是JSON格式）
        content_json <- content(response, "text", encoding = "UTF-8")
        parsed_content <- fromJSON(content_json)
      } else {
        warning(paste("Request failed with status code:", http_status(response)$message))
      }
      tmp <- parsed_content$annotationsByCategory$`pepx-virtual-annotation`
      tmp <- lapply(tmp, function(x){
        x <- paste0(unique(x$cvTermName),collapse = '@')  #以@连接肽段 这些肽段具有重叠的部分
        return(x)
      })
      res[[i]] <- data.frame(Peptide=unlist(tmp),Protein=parsed_content$uniqueName)
    }
    res <- do.call(rbind.data.frame,res)
  }else{
    request_url <- paste0(base_url, "entries/search/peptide.json?peptide=", paste0(Peptide,collapse = ','), 
                          "&no-variant-match=", no_variant_match, 
                          "&mode=", mode)
    # 发送GET请求
    response <- GET(request_url)
    # 检查响应状态码
    if (http_status(response)$category == "Success") {
      print('GET is ok')
      # 解析JSON响应内容（假设API返回的是JSON格式）
      content_json <- content(response, "text", encoding = "UTF-8")
      parsed_content <- fromJSON(content_json)
    } else {
      warning(paste("Request failed with status code:", http_status(response)$message))
    }
    tmp <- parsed_content$annotationsByCategory$`pepx-virtual-annotation`
    tmp <- lapply(tmp, function(x){
      x <- paste0(unique(x$cvTermName),collapse = '@')  #以@连接肽段 这些肽段具有重叠的部分
      return(x)
    })
    res<- data.frame(Peptide=unlist(tmp),Protein=parsed_content$uniqueName)
  }
  return(res)
}



