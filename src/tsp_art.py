from PIL import Image
from PIL import ImageDraw
from stipple import stipple
from tsp import nn,opt2

def draw_route(route, sz, fname, red=set(), green=set(), blue=set()):
    im = Image.new('RGB', sz, (255, 255, 255))
    draw = ImageDraw.Draw(im)
    edges = [[route[i],route[i+1]] for i in range(len(route)-1)] + [(route[-1], route[0])]
    for e in edges:
        draw.line(e, fill=(127, 127, 127), width=2)
    del draw
    im.save(fname)

#以下、「1.画像の読込」に該当
im = Image.open('../img/sample.jpg')
im = im.resize((int(im.size[0]/1.5),int(im.size[1]/1.5)))
gray_im = im.convert('L')
gray_im.save('gray_converted.jpg')
print("横：{}px, 縦:{}px".format(im.size[0],im.size[1]))

#以下、「2. 読込んだ画像を点描画」に該当
points = stipple(gray_im,5)

#以下、「3. 1で生成された点群をTSPとして解く」に該当
nn_route = nn(points)
draw_route(nn_route, im.size, 'tsp_nn.jpg')
route_improved = opt2(nn_route)

#以下、「4. 2で導かれた経路を描画」に該当
draw_route(route_improved, im.size, 'tsp_improved.jpg')

