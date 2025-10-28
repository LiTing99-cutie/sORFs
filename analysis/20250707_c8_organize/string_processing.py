import re

def process_filename(filename):
    """
    处理文件名，去掉数字中的前导零和末尾的.dta扩展名
    
    Args:
        filename (str): 输入的文件名
        
    Returns:
        str: 处理后的文件名
    """
    # 去掉末尾的.dta扩展名
    if filename.endswith('.dta'):
        filename = filename[:-4]
    
    # 使用正则表达式找到数字并去掉前导零
    # 匹配模式：\b0+(\d+)\b 表示单词边界内的前导零+数字
    def remove_leading_zeros(match):
        # 获取匹配的数字部分（去掉前导零）
        number = match.group(1)
        return number
    
    # 应用正则表达式替换
    processed = re.sub(r'\b0+(\d+)\b', remove_leading_zeros, filename)
    
    return processed

# 测试函数
test_filename = "CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.0076.0076.3.dta"
result = process_filename(test_filename)
print(f"原始文件名: {test_filename}")
print(f"处理后: {result}")

# 更多测试用例
test_cases = [
    "CAD20250514licq_BSEP_DDA_60min_21pcw_1_C8_T_T_Slot2-3_1_7020.0076.0076.3.dta",
    "file_001_002_003.dta",
    "sample_0001_0002.dta",
    "test_123_456.dta"
]

print("\n更多测试用例:")
for test in test_cases:
    result = process_filename(test)
    print(f"{test} -> {result}") 