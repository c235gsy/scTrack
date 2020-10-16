from rpy2 import robjects
from rpy2.robjects import numpy2ri
import colour


def get_color_list_R(n, h="c(0, 360) + 15", c=100, l=65):
    r = robjects.r
    r ('''
    library(ggplot2)
    library(scales)
    ''')
    #下面的代码时hue_pal函数的默认设置，
    # 其中h是色相，范围越大，相邻颜色之间差异越大；
    # c是饱和度，值越大色彩越浓艳饱满；
    # l是亮度，大亮小暗
    numpy2ri.activate ()
    comment = "hue_pal(h = {}, c = {}, l = {})({})".format(h, c, l, n)
    color_list = list(r(comment))
    return color_list


def get_color_list_py(n, hue_start=1/24, hue_end=1, saturation=1, luminance=0.65):
    if n <= 0:
        print("n should be bigger than 0")
        return
    elif 0 <= hue_end <= 1 and 0 <= hue_start <= 1:
        out = []
        diss = (hue_end - hue_start) / n
        for i in range(n):
            out.append(colour.Color(hsl=(hue_start + i*diss, saturation, luminance)).hex_l)
        return out


# print(get_color_list_R(3))
# for c in get_color_list_R(3):
#     cc = colour.Color(c)
#     print(cc.hsl, cc.saturation, cc.luminance)
#
#
# print(get_color_list_py(3))
# for c in get_color_list_py(3):
#     cc = colour.Color (c)
#     print (cc.hsl, cc.saturation, cc.luminance)