import requests
url="http://xueshu.baidu.com/usercenter/data/schpaper"
payload={
    "wd":"refpaperuri:(cc56983141c8f43696bd446d45185cb0)",
    "type":"citation",
    "rn":"50",
    "page_no":"1"
}
cookies = dict(cookies_are='working')
r = requests.get(url,params=payload,cookies=cookies)
js=r.json()
resl=js['data']['resultList']
for res in resl:
    authors=res['meta_di_info']['sc_author']
    notus=True
    if authors:
        for author in authors:
            if author['sc_name'][1] =="Zheng-Ming Huang":
                notus=False
        if notus:
            print(res['meta_di_info']['sc_title'])