---
type: regex
target: last_message
pattern: '(\d+(\.\d+)?\s?(mg|mcg|µg|ug|g|ml|mL|IU|iu)\b(\s?/\s?kg)?)|(\d+(\.\d+)?\s?(毫克|毫升|微克|克|片|粒|滴)(\s?/\s?(公斤|千克|kg))?)|((一次|每次|每日|一天)[一二两三四半\d]+\s?(片|粒|袋|支|ml|毫升))'
match: not_contains
weight: 1
---
