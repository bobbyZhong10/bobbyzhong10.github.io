-- 把 pandoc 生成的表格包进 <div class="table-scroll">，以便窄屏横向滚动
function Table(el)
  return pandoc.Div(el, pandoc.Attr("", {"table-scroll"}))
end
