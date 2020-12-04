---
layout: post
title: L0 Gradyan En Küçükleme ile İmge Yumuşatma
slug: l0-gradient-minimization
author: Bahri ABACI
categories:
- Makine Öğrenmesi
- Lineer Cebir
- Nümerik Yöntemler
references: "Image Smoothing via L0 Gradient Minimization"
thumbnail: /assets/post_resources/image_smoothing/thumbnail.png 
---
İmge yumuşatma (Image smoothing) imgeyi oluşturan renk sayısının belirli görüntü işleme teknikleri ile azaltılması işlemidir. Bu işlem özellikle gürültülü görüntüleri gürültüden arındırmakta sıklıkla kullanılan bir yöntemdir. Renklerin azaltılması işlemi bu blogda da değindiğimiz [K-means Kümeleme]({% post_url 2015-08-28-k-means-kumeleme-algoritmasi %}), [Kutu Bulanıklaştırma]({% post_url 2019-02-17-tumlev-imge %}), Gaussian Bulanıklaştırma veya Ortanca Süzgeçler yardımıyla yapılabilir. Bu yazımızda yukarıda sayılan yöntemlerden farklı olarak problemi tamamen matematiksel bir şekilde formülize ederek bir optimizasyon problemine dönüştüren bir yöntem incelenecektir. İncelenecek yöntem, Li Xu ve arkadaşları tarafından 2011 yılında SIGGRAPH konferansında "Image Smoothing via L0 Gradient Minimization" ismi ile yayınlanan bildiriyle literatüre kazandırılmış ve farklı alanlarda kullanılarak binin üzerinde atıf almıştır. 

<!--more-->

Önerilen yöntemin literatürdeki diğer yöntemlerden en önemli farkı; renk yumuşatma işlemini görüntüdeki renk geçilerini sınırlandıracak bir optimizasyon problemi olarak ifade etmesidir. Görüntüdeki her bir renk geçişinin küçük veya büyük bir gradyan yaratacağından önerilen yöntemin temel hedefi görüntüdeki gradyanların sayısının belirli bir limit altında tutulmasıdır. Matematiksel olarak bir değişkenin sayısı $L_0$ norm ile belirlendiğinden, yöntem görüntüdeki gradyanların sayısı $C(S) = \lVert (\partial_x S)^2 + (\partial_y S)^2 \lVert_0$ fonksiyonunu en küçükleyecek $S$ görüntüsünü bulmayı hedeflemektedir. Optimizasyon sonrası bulunacak yumuşatılmış imge $S$ nin, girdi imgesi $I$ ile de benzer olması gerektiğinden incelenen çalışmanın ele aldığı optimizasyon fonksiyonu [Lagrange Çarpanları]({% post_url 2020-01-13-lagrange-carpanlari-yontemi %}) kullanılarak şu şekilde yazılabilir.

$$
E(S) = \lVert S - I \lVert^2 + \lambda C(S)
\tag{1}
\label{objective}
$$

Denklem \ref{objective} de tanımlanan ifadede ilk kısım çıktı imgesi $S$ nin $I$ ile benzerliğini kontrol ederken, ikinci kısım ise çıktı imgesindeki gradyanların sayısını kontrol altında tutmaktadır. Yapılacak fonksiyonun amacı $E(S)$ ile verilen fonksiyonun $S$ üzerinden en küçükleyen $S^\ast$ imgesini bulmaktadır.

Yazılan optimizasyon fonksiyonu oldukça basit görünse de $C(S)$ fonksiyonunun yapısı nedeniyle [Gradyan İniş Yöntemleri]({% post_url 2020-04-08-gradyan-yontemleri-ile-optimizasyon %}) veya diğer gradyan temelli optimizasyon yöntemleri ile doğrudan çözülebilir bir biçimde değildir. Bu nedenle yazarlar 2008 yılında "A new  alternating minimization algorithm for total variation image reconstruction" makalesinde önerilen dönüşümlü optimizasyon stratejisini Denklem \ref{objective} ile verilen probleme uygulayarak sonuca ulşmaya çalışmışlardır. Kullanılan dönüşümlü optimizasyon stratejisinin en önemli özelliği Denklem \ref{objective} ile verilen optimizasyon problemini yardımcı bir alt problem tanımlayarak çözmeye çalışmasıdır. 

Tanımlanan problem $C(S)$ fonksiyonunu $S$ değişkeni yerine sabit kabul edilecek bir yardımcı değişkene bağlı yeniden yazabilmeyi hedeflemektedir. $C(S)$ in, $S$ değişkenine bağlı olmadan yeniden yazılması durumunda Denklem \ref{objective} ifadesi dışbükey olacağından basit bir çözüme sahip olacaktır.

Dönüşümlü optimizasyon stratejisini kullanabilmek için $h=\nabla_x S$ ve $v=\nabla_y S$ tanımlamalarını yaparak Denklem \ref{objective} ile verilen optimizasyon problemini yeniden yazalım.

$$
\begin{aligned}
S^\ast &= \arg \min_{S,h,v} \; E(S,h,v)\\
E(S,h,v) &= \lVert S - I \lVert^2 + \lambda C(h,v) +  \beta \left( \lVert \nabla_x S - h \lVert^2 + \lVert \nabla_y S - v \lVert^2 \right)
\end{aligned}
\tag{2}
\label{objective_aux}
$$

Yazılan ifadede $C(h,v) = \lVert h^2 + v^2 \lVert_0$ gradyan sayma fonksiyonunu göstermektedir ve $\beta$ kullanılan yardımcı değişkenlerin orjinal değişkenlere yakın olmasını zorlamaktadır. Denklem \ref{objective_aux} incelendiğinde $\beta$ katsayısının çok büyük seçilmesi durumunda $h=\nabla_x S$ ve $v=\nabla_y S$ şartı sağlanacak ve Denklem \ref{objective} ile verilen orjinal probleme yakınsayacaktır. 

Denklem \ref{objective_aux} ile yazılan optimizasyon probleminin Denklem \ref{objective} den en önemli farkı, $C$ fonskiyonunu $S$ değişkeninden bağımsız bir şekilde yeniden yazması ve $S$ üzerinden optimizizasyonu oldukça kolaylaştırmasıdır. Ancak bunu yaparken denkleme eklediği iki yeni değişken $(h,v)$ üzerinden de yeni bir optimizasyon problemi ortaya çıkmaktadır. Denklem \ref{objective_aux} ile verilen problemin çözümü için $(h,v)$ üzerinden optimizasyon problemiyle $S$ üzerinden optimizasyon problemi dönüşümlü olarak çözülmelidir.

### Alt Problem 1: $h,v$ bilinmesi durumunda $S^\ast$

Denklem \ref{objective_aux} ile verilen $E(S,h,v)$ ifadesinin $S$ e bağlı en küçük değeri, $S$ değişkenine göre türev alınarak bulunabilir. Bu ifadede $h,v$ değerleri $S$ değerine bağlı olmadığından türev işlemi sırasında skalar olarak ele alınacağından, aşağıdaki çözüm elde edilir.

$$
\begin{aligned}
\frac{\partial E(S,h,v)}{\partial S} &= 2(S - I) + 2 \beta \left( \nabla^\intercal_x(\nabla_x S - h) + \nabla^\intercal_y(\nabla_y S - v) \right) = 0\\
\Rightarrow & I + \beta \left( \nabla^\intercal_x h + \nabla^\intercal_y v \right) = S + \beta \left( \nabla^2_x S + \nabla^2_y S \right)\\
\Rightarrow & I + \beta \left( \nabla^\intercal_x h + \nabla^\intercal_y v \right) = S \left(1 + \beta \left( \nabla^2_x + \nabla^2_y \right) \right)\\
\Rightarrow & B = S \left(1 + \beta \Delta \right)
\end{aligned}
$$

Elde edilen ifadede sol tarafta $B=I + \beta \left( \nabla^\intercal_x h + \nabla^\intercal_y v \right)$ ile gösterilen değişkende bilinen değerler toplanmış, sağ tarafta ise bilinmeyen $S$ değeri bırakılmıştır. Denklemin sağında yer alan $\Delta = \nabla^2_x + \nabla^2_y$ ifadesi Laplace operatörüdür ve [Poisson Denklemi Yardımıyla Görüntü Düzenleme]({% post_url 2020-09-20-poisson-denklemi-yardimiyla-goruntu-duzenleme %}) yazımızda da değindiğimiz üzere iki boyutlu imgeler için 

$$\nabla^2_x + \nabla^2_y = \begin{bmatrix}\phantom{+}0 & -1 & \phantom{+}0\\-1 &\phantom{+}4 & -1\\\phantom{+}0 & -1 & \phantom{+}0\end{bmatrix}$$ 

çekirdeği ile evrişim işlemi sonucunda hesaplanır. Benzer şekilde $\nabla_x = \begin{bmatrix}-1 & 1\end{bmatrix}$ ve $\nabla_y = \begin{bmatrix}-1 && 1\end{bmatrix}^\intercal$ ifadeleri de yönlü gradyan operatörünü göstermektedir. 

[Poisson Denklemi Yardımıyla Görüntü Düzenleme]({% post_url 2020-09-20-poisson-denklemi-yardimiyla-goruntu-duzenleme %}) yazımızda yaptığımız şekilde bu denklemleri bir matris içerisinde yazmak istersek aşağıdaki matris ifadesi elde edilir.

$$
\begin{bmatrix}
\dots \\ \dots \\ \dots \\
\hline
\dots \\ B(x,y) \\ \dots \\ 
\hline
\dots \\ \dots \\ \dots
\end{bmatrix}=
\underbrace{
\begin{bmatrix}
&&&&\dots\\
&&&&\dots\\
&&&&\dots\\
\hline
&&&&\dots\\
0 \dots 0 & -\beta & 0 \dots 0 & -\beta & 1 + 4\beta & -\beta & 0 \dots 0 & -\beta & 0 \dots 0\\
&&&&\dots\\
\hline
&&&&\dots \\
&&&&\dots \\
&&&&\dots \\
\end{bmatrix}
}_{\mathbf{P}}
\begin{bmatrix}
\dots \\ S(x,y-1) \\ \dots\\
\hline
S(x-1,y) \\ S(x,y) \\ S(x+1,y)\\ 
\hline
\dots \\S(x,y+1) \\ \dots
\end{bmatrix}
\tag{3}
\label{s_solution}
$$

Elde edilen $B=PS$ ifadesinde $P$ Poisson matrisi olarak bilinmektedir. Denklem doğrusal bir denklem takımı olduğundan imgedeki tüm pikseller için $B=PS$ denklemleri yazılarak $S^\ast=B^{-1} P$ işlemi ile $S^\ast$ imgesi çözülebilir. $P$ matrisinin oldukça seyrek bir matris olması özelliği göz önüne alındığında işlemler  [Successive Over-Relaxation]({% post_url 2020-09-20-poisson-denklemi-yardimiyla-goruntu-duzenleme %}) yöntemi kullanılarak oldukça hızlı bir şekilde çözülebilir.


### Alt Problem 2: $S^\ast$ bilinmesi durumunda $h,v$

Yazımızın ilk kısmında da değindiğimiz üzere Denklem \ref{objective_aux} ile yazılan problem $S$ üzerinden en küçükleme işlemini oldukça kolaylaştırmaktadır. Ancak $S^\ast$ sabit kabul edilmesi durumunda yazılacak olan

$$
E(S^\ast,h,v) = \lambda C(h,v) + \beta \left( \lVert \nabla_x S - h \lVert^2 + \lVert \nabla_y S - v \lVert^2 \right)
$$

ifadesinin çözümü dış bükey veya türevlenebilir olmadığından problemin nümerik olarak çözülmesi hala problemlidir. Ancak basit bir akıl yürütme yapılarak $E(S^\ast,h,v)$ ifadesini en küçükleyen $h,v$ değerleri bulunabilir. Bunun için öncelikle matris formunda verilen ifadeyi $p=(x,y)$ gibi seçilen herhangi bir piksel için yazıp, yazılan tüm ifadeyi $\beta > 0$ olmak üzere $\beta$ ile normalize edersek

$$
E_p(S^\ast,h_p,v_p) = \frac{\lambda}{\beta} C(h_p,v_p) +  \left( \lVert \nabla_x S_p - h_p \lVert^2 + \lVert \nabla_y S_p - v_p \lVert^2 \right)
\tag{4}
\label{hv_solution}
$$

Eşitliği elde edilir. Denklem \ref{hv_solution} te tanımlanan ifadede $C(h_p,v_p)$ fonksiyonu da parçalı fonksiyon kullanılarak aşağıdaki şekilde yazılabilir.

$$
C(h_p,v_p) = 
\begin{cases}
1, & h_p \neq 0 \text{ veya } v_p  \neq 0\\
0, & h_p = 0 \text{ ve } v_p  = 0
\end{cases}
$$

Denklem \ref{hv_solution} ile verilen ifade ve $C$ fonksiyonu birlikte ele alındığında $h_p=0$ ve $v_p=0$ seçilmesi durumunda $E(S^\ast, 0,0) = (\nabla_x S_p)^2 + (\nabla_y S_p)^2$ şeklinde hesaplanacaktır. $h_p,v_p$ değerlerinin sıfırdan farklı seçilmesi durumunda da en mantıklı seçim ikinci terimleri sıfır yapacak olan $h_p=\nabla_x S_p$ ve $v_p=\nabla_y S_p$ seçilmesi olacaktır. Bu drumumda ise toplam hata $E(S^\ast, \nabla_x S_p,\nabla_y S_p) = \frac{lambda}{\beta}$ şeklinde bulunacaktır.

Bu durumda $(\nabla_x S_p)^2 + (\nabla_y S_p)^2$ gradyan toplamı $\frac{lambda}{\beta}$ değerinden büyük olması durumunda $h_p=\nabla_x S_p, \; v_p=\nabla_y S_p$ seçilmesi mantıklıyken, diğer durumda $h_p=0, \; v_p=0$ en iyi seçim olacaktır. Bu seçim stratejisi aşağıdaki parçalı fonskiyon ile gösterilebilir.

$$
(h_p,v_p) = 
\begin{cases}
(0,0), &  (\nabla_x S_p)^2 + (\nabla_y S_p)^2 \leq \frac{lambda}{\beta} \\
(\nabla_x S_p, \nabla_y S_p), & \text{ diğer}
\end{cases}
\tag{5}
\label{hp_solution_partial}
$$

Yukarıda elde edilen iki alt problem ve çözümü birlikte ele alındığında Denklem \ref{objective_aux} ile tanımlanan problemin çözümü aşağıdaki algoritma yardımıyla hesaplanır.

> - **GİRDİLER**
>   - $I$: imge
>   - $\beta, \beta_0, \beta_\text{max}$: yumuşatma parametresi
>   - $\lambda$: gradyan cezası
>   - $\kappa$: $\beta$ adım büyüklüğü 
> - **ÇIKTILAR**
>   - $S$: yumuşatılmış imge
>
> *****************************
>
> - $S$ = $I$
> - $\beta = \beta_0$
> - **while** $\beta < \beta_\text{max}$
>   - Denklem \ref{hp_solution_partial} ve $S$ yi kullanarak $h,v$ değerlerini hesapla
>   - Bulunan $h,v$ değerlerini kullanarak $B=I + \beta \left( \nabla^\intercal_x h + \nabla^\intercal_y v \right)$ matrisini hesapla
>   - $B$ ve Denklem \ref{s_solution} ü kullanarak yeni $S$ yi hesapla
>   - $\beta = \kappa \beta$ ile $\beta$ değerini artır
> - **return** $S$

Verilen algoritmadan görüldüğü üzere ilk olarak girdi imgesi $S$ imgesine atanarak ilklendirme yapılmakta ve bu $S$ değeri kullanılarak $h,v$ gradyanları bulunmakta. Ardından bu değerler kullanılarak bilinenler matrisi $B$ hesaplanmaktadır. $B$ bulunduktan sonra Denklem \ref{s_solution} ile yazılan Poisson eşitliği çözülerek $S$ imgesi bulunmaktadır. Her iterasyonda bir önceki iterasyonda kullanılan $\beta$ değeri $\kappa=2.0$ gibi sabit bir sayı ile çarpılarak ertırılmakta ve $\beta$ değerinin olması gerektiği gibi çok yüksek bir sayıya ulaşması sağlanmaktadır.

Yukarıda verilen algoritma IMLAB görüntü işleme kütüphanesi kullanılarak `L0Minimize(matrix_t *input, float lambda, float beta0, float betaMax, float kappa, matrix_t *output)` fonskiyonu şeklinde yazılmıştır. Yazılan fonksiyon makalede de önerilen parametreler kullanılarak aşağıdaki şekilde kullanılmıştır.

```c
matrix_t *img = imread("..//data//flower.bmp");
matrix_t *output = matrix_create(uint8_t);

// set the L0 minimization paramaters
float lambda = 0.03f;
float kappa = 2.0f;
float beta0 = 2.0f * lambda;
float betaMax = 1e5f;

// create L0 regularized image
L0Minimize(img, lambda, beta0, betaMax, kappa, output);
```

Yazılan kod parçası farklı imgeler üzerinde çalıştırılarak aşağıdaki sonuçlar elde edilmiştir. Üretilen sonuçların giriş kısmında değinilen benzer algoritmalarla karşılatırılabilmesi için Gaussian Bulanıklaştırma $(11\times 11)$, Ortanca Süzgeç $(11\times 11)$ ve Seçici Gaussian Bulanıklaştırma $(11\times 11, s=0.4)$ gibi farklı yöntemlerin çıktıları da verilmiştir.

| Kaynak İmge | Gaussian Bulanıklaştırma | Ortanca Süzgeç | Seçici Gaussian Bulanıklaştırma | L$0$ Yumuşatma
:-------:|:----:|:----:|:---:|:---:|
![Image Smoothing][pattern] | ![Gaussian Bulanıklaştırma][pattern_gaussian] | ![Ortanca Süzgeç][pattern_median] | ![Seçici Gaussian Süzgeç][pattern_selective] | ![L0 Image Smoothing][pattern_result] |
![Image Smoothing][flower] | ![Gaussian Bulanıklaştırma][flower_gaussian] | ![Ortanca Süzgeç][flower_median] | ![Seçici Gaussian Süzgeç][flower_selective] | ![L0 Image Smoothing][flower_result] |
![Image Smoothing][lanterns] | ![Gaussian Bulanıklaştırma][lanterns_gaussian] | ![Ortanca Süzgeç][lanterns_median] | ![Seçici Gaussian Süzgeç][lanterns_selective] | ![L0 Image Smoothing][lanterns_result] |

Verilen tablo incelendiğinde önerilen yöntemin pencere boyundan bağımsız olarak çok geniş alanlarda da yumuşatma işlemi yapabildiği ve bunu yaparken görüntüyü herhangi bir şekilde bulanıklaştırmadığı görülmektedir. 

Yazıda yer alan analizlerin yapıldığı kod parçaları, görseller ve kullanılan veri setlerine [Image Smoothing](https://github.com/cescript/image_smoothing) GitHub sayfası üzerinden erişebilirsiniz.

**Referanslar**
* Xu, Li, et al. "Image smoothing via L 0 gradient minimization." Proceedings of the 2011 SIGGRAPH Asia Conference. 2011.
* Wang, Yilun, et al. "A new alternating minimization algorithm for total variation image reconstruction." SIAM Journal on Imaging Sciences 1.3 (2008): 248-272.

[RESOURCES]: # (List of the resources used by the blog post)
[pattern]: /assets/post_resources/image_smoothing/pattern.png
[pattern_gaussian]: /assets/post_resources/image_smoothing/pattern_gaussian.png
[pattern_median]: /assets/post_resources/image_smoothing/pattern_median.png
[pattern_selective]: /assets/post_resources/image_smoothing/pattern_selective.png
[pattern_result]: /assets/post_resources/image_smoothing/pattern_result.png

[flower]: /assets/post_resources/image_smoothing/flower.png
[flower_gaussian]: /assets/post_resources/image_smoothing/flower_gaussian.png
[flower_median]: /assets/post_resources/image_smoothing/flower_median.png
[flower_selective]: /assets/post_resources/image_smoothing/flower_selective.png
[flower_result]: /assets/post_resources/image_smoothing/flower_result.png

[lanterns]: /assets/post_resources/image_smoothing/lanterns.png
[lanterns_gaussian]: /assets/post_resources/image_smoothing/lanterns_gaussian.png
[lanterns_median]: /assets/post_resources/image_smoothing/lanterns_median.png
[lanterns_selective]: /assets/post_resources/image_smoothing/lanterns_selective.png
[lanterns_result]: /assets/post_resources/image_smoothing/lanterns_result.png