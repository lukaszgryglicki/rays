# rays
NURBS and other objects raytracing, written in C

# descr

Only in Polish:




Politechnika Warszawska
Wydział Matematyki i Nauk Informacyjnych
kierunek: Informatyka
specjalizacja: CAD/CAM





Łukasz Gryglicki
Numer indeksu: 181556




Praca magisterska

Raytracing powierzchni NURBS






Promotor pracy magisterskiej
prof. Krzysztof Marciniak


Mińsk Mazowiecki, kwiecień 2006

Spis treści
	1.Wstęp
		1.1. Uproszczony algorytm raytracingu [1]
	2. Powierzchnie NURBS
		2.1 Definicja powierzchni NURBS
			2.1.1	Algorytm deCastlejau:
			2.1.2	Obliczenie NURBS dla zadanych u,v:
		2.2  Problemy dotyczące raytracingu powierzchni 			NURBS
		2.3 Metoda tradycyjna raytracingu powierzchni 				NURBS
		2.4 Metoda raytracingu pow. NURBS przez 					triangulację
			2.4.1 Obliczanie normalnych powierzchni NURBS
			2.4.2 Triangulacja, sposób podziału powierzchni
			2.4.3 Struktura trójkąta GPT
			2.4.4 Obliczanie powierzchni w której zawiera się trójkąt
			2.4.5 Obliczanie normalnych trójkąta na podstawie 						wierzchołków
			2.4.6 Algorytmy interpolacji
			2.4.7 Algorytm przecięcia linii
		2.5 Wady obliczania raytracingu przez triangulację
		2.6 Zalety przyjętego przeze mnie rozwiązania
	3. Główne algorytmy
		3.1 Algorytmy przecięć
		3.2 Hierarchia  obiektów otaczających (BVH)
			3.2.1 Algorytm generacji drzewa AABB
				3.2.1.1 Algorytm generujący drzewo z listy
				3.2.1.2 Algorytm minimalizacji pełnej
				3.2.1.3 Algorytm minimalizacji częściowej
				3.2.1.4 Algorytm szybki
				3.2.1.5 Algorytm wokselowy
				3.2.1.6 Algorytm sasiedztwa
			3.2.2 Zapis i odczyt drzewa AABB
				3.2.2.1 Zapis drzewa AABB
				3.2.2.2 Odczyt drzewa AABB
			3.2.3 Algorytm przecięcia (Intersection)
		3.3 Algorytm anty-aliasingu i rzucania cieni
		3.4 Efekt Fresnela
		3.5 Algorytm wyznaczania koloru tła
		3.6 Algorytm raytracingu
			3.6.1 Główna funkcja raytracingu – kolor rekurencyjny
	4. Opis struktur danych oraz ich zastosowania
		4.1 Struktury danych
		4.2 Zastosowanie struktur danych w algorytmach
	5. Zaimplementowane funkcjonalności, programy
		5.1 Odczyt i zapis danych
			5.1.1 Odczyt danych
				5.1.1.1 Wczytywanie definicji sceny
				5.1.1.2 Wczytywanie tekstur i scen wygenerowanych cczęśściowo
				5.1.1.3 Wczytywanie drzewa AABB
			5.1.2. Zapis danych
				5.1.2.1 Zapis obrazów wynikowych/częściowych
				5.1.2.2 Zapis binarnych scen
				5.1.2.3 Zapis obliczonego drzewa AABB
				5.1.2.4 Zapis przekształceń w trybie podglądu
				5.1.2.5 Konwersje międzyformatowe
		5.2 Przekształcenia sceny i kopiowanie obiektów
			5.2.1 Przekształcenia obiektów
				5.2.1.1 Przekształcenia pojedyńczego trójkąta
				5.2.1.2 Przekształcenia list trójkątów.
				5.2.1.3 Przekształcenia list powierzchni NURBS
			5.2.2 Kopiowanie obiektów
				5.2.2.1 Kopiowanie pojedynczego trójkąta
				5.2.2.2 Kopiowanie listy trójkątów
				5.2.2.3 Kopiowanie powierzchni NURBS
		5.3 Inne algorytmy
			5.3.1 Obsługa komentarzy w pliku wejściowym
			5.3.2 Obsługa sterowania przez internet
			5.3.3 OpenGL GUI
				5.3.3.1 Wyświetlanie raytracingu w OpenGL
				5.3.3.2 Podgląd sceny w OpenGL
			5.3.4 Chropowatość powierzchni
			5.3.5 Obsługa sygnału błedu segmentacji
			5.3.6. Generowanie sceny w odcieniach szarości
			5.3.7 Łączenie i transformacje AABB-drzew
				5.3.7.1 Algorytm łączenia AABB-drzew:
				5.3.7.2 Algorytm przekształcania AABB-drzew
	6.Opis użycia wszystkich programów, dane techniczne
		6.1 Dane techniczne
		6.2 Opisy programów
			6.2.1 Program RAYS
				6.2.1.1 Opis opcji wiersza poleceń
				6.2.1.2 Opis poleceń w oknie graficznym OpenGL
				6.2.1.3 Opis poleceń w oknie podglądu OpenGL
				6.2.1.4 Opis poleceń serwera rays.
			6.2.2 Program NURBS
				6.2.2.1 Algorytm interpolacji globalnej
				6.2.2.2 Opcje programu NURBS
				6.2.2.3 Używanie programu NURBS
			6.2.3 Program NURBS2DAT
			6.2.4 Program IGES2DAT
			6.2.5 Program ULI2DAT
			6.2.6 Program MD22DAT
			6.2.7 Program 3DS2TRI
			6.2.8 Program 3DS2DAT
			6.2.9 Program 60FACES
			6.2.10 Program TORRUSGEN
			6.2.11 Program CONE
			6.2.12 Program CUBE
			6.2.13 Program CYLINDER
			6.2.14 Program BALL
			6.2.15 Program RANDNURB/RANDNURBFULL
			6.2.16 Program RTRIANGLE
			6.2.17 Program TABLE
			6.2.18 Program TERMINAL
			6.2.19 Program GETBMP
			6.2.20 Program TEX
			6.2.21 Program WRAPPER
			6.2.22 Program BTREECONV
	7. Przykładowe sceny, sposób renderowania, metody kompilacji
		7.1 Przykładowe sceny
		7.2 Kilka przykładów renderowania
		7.3 Sposoby kompilacji
		7.4 Struktura katalogów
		7.5 Przykładowe obrazy
	8. Bibliografia












	1.	Wstęp

	Raytracing (ang. śledzenie promieni) jest techniką renderowania realistycznych scen 3D.W przeciwieństwie do standardowych metod wyświetlania sceny w grafice komputerowej, które dokonują szeregu uproszczeń, raytracing oblicza dokładny wygląd sceny poprzez prześledzenie biegu każdego promienia od oka obserwatora do wybranych pikseli ekranu. Uwzględniane są przy tym zjawiska fizyczne takie jak: efekt Fresnela, ugięcie światła, załamanie, pochłonięcie przez materiał. W przeciwieństwie do algorytmów uproszczonych (używanych w grafice komputerowej, np. w grach) raytracing na ogól nie daje się obliczyć w czasie rzeczywistym na współczesnym sprzęcie.

	Są próby implementacji raytracingu bezpośrednio na procesorach graficznych, używając do tego dodatkowego przejścia renderowania: pixel/vertex shader'ów. Współczesne karty graficzne implementujące Pixel Shader 3.0 (lub conajmniej 2.0) są w stanie generować kilkanaście klatek na sekundę, przy małej rozdzielczości (400x300, 640x480) i ograniczonej rekursji (do co najwyżej kilku odbić/załamań), takie wyniki są osiągane przez najszybsze/najlepsze programy tego typu dostępne na rynku.

  	Przykładowy system http://gpurt.sourceforge.net/DA07_0405_Ray_Tracing_on_GPU-1.0.5.pdf, wymaga karty graficznej GeForce 6600 GT (PS 3.0), osiąga wynik 0.25 FPS dla średniej wielkości scen. Jednak możliwości programowania jednostek PS i VS są jeszcze zbyt małe aby była możliwość zaimplementowania tam odpowiednich procedur przecinających powierzchnię Spline/NURBS z promieniem, poza tym zwolniłoby to rendering o następne kilka rzędów wielkości.

	Tworząc algorytmy raytracingu zrezygnowałem więc z działania w czasie rzeczywistym, skoncentrowałem się na dokładnym (nie koniecznie najszybszym) zaimplementowaniu wszystkich efektów raytracingu. Zastosowałem także szczególne podejście do powierzchni NURBS, które opiszę poniżej. Mój algorytm daje dobre efekty, istnieją jednak lepsze algorytmu, generujące na przykład “miękkie cienie” oraz likwidujące efekty aliasingu, czy pozwalające na rozszczepienie światła.


	Takim algorytmem jest na przykład algorytm : “Forward Raytracing” [3] - jest to metoda śledzenia promieni, która w przeciwieństwie do normalnego (czyli wstecznego – backward) raytracingu  zaczyna się nie od kamery, lecz od źródeł światła. Określenie “forward” oznacza, że symulowane promienie poruszają się w tym samym kierunku co rzeczywistym świecie. Daje ona lepsze wyniki niż raytracing wsteczny, lecz wymaga nieporównywalnie więcej mocy obliczeniowej, więc jest używana tylko do pewnych specjalnych zastosowań, mających niewiele wspólnego z fotorealistycznym renderingiem, takich jak badanie właściwości sprzętu optycznego oraz w rozwiązaniach hybrydowych w połączeniu ze zwykłym raytracingiem. Metoda wymaga dużo większej mocy obliczeniowej, ponieważ potrzebujemy wystrzelić promienie od źródła światła we wszystkich kierunkach, a tylko niewielka część z nich dotrze do oka obserwatora (szansa jest bardzo mała), poza tym musimy zrobić to dla każdego źródła światła. Zaletą jest to, że promienie poruszają się tak w rzeczywistości co pozwala uzyskać pewne efekty nieosiągalne zwykłym raytracingiem, np. efekty soczewkowe: (rozszczepienie światła, ogniska świateł skupionych przez soczewkę, miękkie cienie, światła odbite - “zajączki”). Poza tym w raytracingu wstecznym (backward), określanie czy w danym punkcie jest cień odbywa się na podstawie testu czy z danego punktu wysyłając promień do źródła światła, przetnie on jakiś obiekt czy nie. Powoduje to, że granice cienia są “ostre” a cienie “zupełnie czarne”. Innym problemem raytracingu (a także innych metod generowania obrazów) jest aliasing (ostre krawędzie, artefakty na teksturach itp). W moim algorytmie jest możliwość ustawienia prostego, programowego (software) antyaliasingu 2x2. Spowalnia to proces około 4x.

	Inne metody generowania obrazów to: radiosity, photon map, ray-casting [1] (będący raytracingiem o rekursji ograniczonej do 0 oraz różne kombinacje wyżej wymienionych.

	W moim algorytmie próbuję poradzić sobie z różnymi problemami raytracingu: takimi jak: rozszczepienie światła, aliasing, szybkość, działania, teksturowanie, cienie (problem czarnego cienia i problem zbyt ostrego cienia), chropowate powierzchnie, efekt Fresnela, głębokość rekursji.

	Raytracing, mimo swoich zalet, nie jest idealnym sposobem tworzenia obrazów. Przede wszystkim słabo radzi sobie ze światłem rozproszonym i z modelowaniem bardziej skomplikowanych źródeł światła (w moim modelu zakładam, że źródła światła są punktami, co nie zdarza się w praktyce). Efektem są bardzo ostre, nierealistycznie wyglądające, krawędzie cieni.

	Ponieważ raytracing operuje na pojedynczych promieniach, nie możne prawidłowo modelować dyfrakcji, np. Na pryzmacie, interferencji fal świetlnych i innych zjawisk falowych. Pewnym rozwiązaniem jest rozdzielenie pojedynczego promienia na kilka promieni reprezentujących różne długości fal (w rzeczywistości jest to problem ciągły a nie dyskretny). Wybrałem podejście jak najbardziej zbliżone do modelu maszynowego. Kolor jest zapisywany w obrazach jako RGB lub RGBA. Model RGB (24bity) definiuje po 8 bitów na czerwony, zielony i niebieski kolor. Odpowiednie kombinacje tych składowych dają każdy widzialny kolor. Założyłem, że przez jeden pixel wystrzelane są 3 promienie: czerwony, zielony i niebieski. Każdy z nich może mieć natężenie 8bitowe, tj: 0-255.

	Stworzono wiele technik które stosowane samodzielnie, bądź razem z raytracingiem pozwalają obejść te wady. Są to m.in. radiosity, photon map, global ilumination [5], forward raytracing.

1.1. Uproszczony algorytm raytracingu [1]

      Z punktu, w którym znajduje się kamera wypuszczany jest promień (półprosta) w kierunku rzutni. Rzutnia podzielona jest na piksele, jeden (lub więcej) promieni przechodzi przez każdy pixel Typowe rzutnie mają 800x600 pikseli, co daje 480000 promieni, poza tym oddzielne promienie są używane dla kolorów R,G,B: daje to razem 1440000 promieni.

	1.	Wyszukiwane są wszystkie przecięcia promienia z obiektami. Tutaj szybkość algorytmu zależy nie tyle od efektywności wyznaczania współrzędnych punktu przecięcia promienia z powierzchnią obiektu, ale od tego by tych dokładnych obliczeń było jak najmniej, by przetwarzać tylko te obiekty, które wysłany w p.1 promień przecina. Mówiąc obrazowo: jeśli scena zawiera milion obiektów, by za każdym razem nie testować milion razy, czy promień przecina każdy z obiektów, ale testować jedynie te, które potencjalnie mogą mieć punkt przecięcia z promieniem. Ważne więc jest zastosowanie odpowiedniej lokalizacji (wykrywanie kolizji). Należy utworzyć drzewo hierarchiczne (BV – Bounding Volume)

	1.	Spośród uzyskanych punktów przecięć wybiera się ten, który leży najbliżej kamery.

	1.	Punkt ten jest następnie przetwarzany. Najpierw wypuszczamy promień z tego punktu w kierunku światła na scenie, by określić czy oświetla przetwarzany punkt. Na tym etapie można wyznaczyć cienie, testując czy odcinek pomiędzy punktem przecięcia, a światłem przecina jakiś obiekt - innymi słowy, czy jakiś obiekt zasłania światło. Tutaj możemy też uwzględnić jak “głęboki będzie ten cień”, na podstawie współczynnika przezroczystości obiektów na drodze do źródła światła. Następnie oblicza się, używając zadanego modelu oświetlenia (np. Lamberta, Phonga, Metal Shading), jasność punktu. Dodatkowo uwzględnia się takie parametry jak kolor punktu: np. Oczytany z tekstury lub zadany jako własność materiału

	1.	Jeśli obiekt jest przezroczysty lub odbija światło, to z tego punktu wypuszczamy dodatkowe promienie może to być zarówno promień odbity, jak i promień załamany - dla tych promieni algorytm jest powtarzany od punktu 2. Wówczas, nim przypisze się kolor danemu pikselowi, przetwarzane jest drzewo promieni. Maksymalną głębokość rekursji można ograniczyć, gdyż mogłaby być nieskończona, np. wewnątrz sześcianu zbudowanego z luster.

	1.	Powierzchnie NURBS
	1	 Definicja powierzchni NURBS
	Powierzchnię NURBS (Non-Uniform Rational Bspline) możemy zdefiniować następująco:
	Wprowadźmy u i v - zmienne niezależne z przedziału [0,1].
u i v będą parametryzacją powierzchni NURBS: [0,1] x [0,1]. Powierzchnię w punkcie u,v będziemy oznaczać P(u,v). u,v e [0,1]
	Wprowadźmy p1 i p2 – wymiary powierzchni NURBS (odpowiednio w kierunku u i v), p1,p2 >=1
	Wprowadźmy n1 i n2: ilość punktów kontrolnych w kierunku u i v. n1> p1, n2>p2. Punkty kontrolne tworzą macierz M[n1 x n2]
	Wprowadźmy węzły (nodes): t1[n1], t2[t2]. Podział wartości w węzłach jest najczęściej jednorodny, np. dla n1=5 mamy t1=[0 0.25 0.5 0.75 1]. Użytkownik może oczywiście wprowadzić własne węzły (nodes)
	Wprowadźmy: m1, m2: m1=n1+p1, m2=n2+p2
	Wprowadźmy węzły (knots): knot1[m1+1], knot2[m2+1], np.: dla wymiaru = 3 i ilości punktów = 4 mamy p1=3, n1=4, m1=7 knot1[0 0 0 0 1 1 1 1] – typowe wartości węzłów. Węzły (knots) mają większy wymiar aby rekurencyjna funkcja bazowa bspline mogła “sięgać” do elementów przesuniętych o wartości nie przekraczające wymiaru (p). Jest wiele metod na obliczenie węzłów (knots), najczęściej oblicza się je używając węzłów (nodes). Powszechne metody parametryzacji to: jednorodna, naturalna, długość łuku [7] p+1 pierwszych elementów knot to 0, p+1 ostatnich to 1, a pozostałe to rozkład równomierny np. dla p=2 i n = 5 mamy: t[0 0.5 1], knot[0 0 0 0.25 0.5 0.75 1 1 1]

Teraz możemy wprowadzić funkcję bazową bspline:

	Rekurencyjny wzór na obliczenie funkcji bazowej jest następujący [11]:





	Wtedy punkty krzywej NURBS wyrażają się następującym wzorem [11]:

W moim algorytmie używam:
	b0(knot, t, i): jeżeli  t >= knot[i]  oraz t <= knot[i+1] zwróć 1, inaczej zwróć 0, gdzie t to wartość parametru z [0,1]a knot to wprowadzone powyżej węzły. Funkcja ta jest używana na najniższym poziomie rekursji (stopnia zerowego). Aby obliczyć dla p-tego stopnia należy zbudować drzewo deCastlejau. Innym rozwiązaniem jest obliczenie rekurencyjne, w którym dla p-tego wymiaru wywołujemy funkcję bazową 2^p razy (ale jest to nieoptymalne). W mojej pracy wykorzystuje algorytm deCastlejau o złożoności O(N^2).

	2.1.1	Algorytm deCastlejau:
Drzewo(knot, t, p, n):
Wejście: knot[], t, p, n: węzły, wartość t, wymiar i ilość punktów
Wyjście: tablica (w której jest pełne drzewo deCastlejau)
od i = 0 do n wykonaj:
	od j = 0 do p wykonaj:
		jeśli j = 0
			od k = i do i + p wykonaj: tablica[i][j][k-i] = b0(knot, t, k)
		wpp:
			od k = i do i + p wykonaj
			          up = (t - knot[k]) * tablica[i][j-1][k-i];
		      	          down = knot[k+j] – knot[k];
			          jezeli down <> 0 f1 = up/down wpp f1 = 0
		   	          up = (knot[k+j+1] - t) * tablica[i][j-1][k-i+1];
		                        down = knot[k+j+1] – knot[k+1];
			          jezeli down <> 0 f2 = up/down wpp f2 = 0
	  		          tablica[i][j][k-i] = f1 + f2;

Po wykonaniu tego algorytmu mamy pełne drzewo deCastlejau dla każdego wymiaru <= p, dla punktów 0..n

	2.1.2	Obliczenie NURBS dla zadanych u,v:

NURBS(u,v)

Wejście: u,v
Wyjście: wartość powierzchni NURBS w zadanych u,v


D[i][j] – punkty kontrolne powierzchni NURBS
w[i][j] – wagi punktów: 0 (nie brany pod uwagę) im większa tym bliżej tego punktu zostanie “dociągnięta” powierzchnia NURBS. Gdy wszystkie wagi = 1, to powierzchnia BSpline

Obliczenie drzew deCastlejau dla knot1 (kierunek U) i knot2 (kierunek V)

 tab1 = Drzewo(knot1, u, p1, n1)
 tab2 = Drzewo(knot2, v, p2, n2)

 down = up = 0;
od i = 0 do n1 wykonaj:
	tutaj mamy wyliczenie dla wymiaru p1, korzeń drzewa w elemencie o indeksie 0
	b1 = tab1[i][p1][0];
	dla j = 0 do n2 wykonaj:
		b2 = tab2[j][p2][0];
        		factor = D[i][j]*b1*b2*w[i][j];
	              up += factor;
		down += b1*b2*w[i][j];
 jeżeli down <> 0 zwróć up/down wpp. zwróć 0


	1	 Problemy dotyczące raytracingu powierzchni NURBS

	Bezpośredni raytracing powierzchni NURBS jest trudnym problemem. Są dwa główne powody dla których unika się bezpośredniego przetwarzania powierzchni Spline/NURBS:

	a) Szybkość działania. Metody obliczania punktu przecięcia powierzchni NURBS i promienia są iteracyjne, ich zbieżność jest bardzo wolna, często są problemy z tą zbieżnością, czas działania jest o kilka rzędów większy niż przy działaniu na trójkątach. Z testów wynika, że w najlepszym przypadku (dosyć płaska powierzchnia NURBS, ustawiona pod kątek normalnym do obserwatora) jest co najmniej 30x wolniej (Rzadko osiąga się to 30x, częściej jest to 100x). W sytuacji typowej jest to 100x-1000x wolniej!


	b) Niestabilność. Ogólnie są bardzo duże problemy ze zbieżnością metod numerycznych, wiele bibliotek/programów implementujących raytracing NURBS ma z tym problemy. Jest dużo przypadków szczególnych, należy robić wstępne podziały powierzchni na płatki wypukłe metoda deCastlejau, poza tym te płatki muszą mieć odpowiednio mała krzywizną. Promień musi być reprezentowany jako dwie przecinające się płaszczyzny (promień jest krawędzią ich przecięcia). Prawie żaden program NIE działa z dowolnymi powierzchniami NURBS: samoprzecięcia powierzchni, wielokrotne przecięcia z powierzchnią, rozbieżność Newtona, NURBS'y mogą mieć kawałki z ciągłością C0, w takich punktach obliczenia numeryczne zawodzą.

Zacytuje: [4]  “Using the direct raytracing NURBS surfaces, one can achieve better quality of rendered images. Although, many such approaches have already been presented, almost all of them suffer from numerical problems or do not work in some special cases. “

2.3 Metoda tradycyjna raytracingu powierzchni NURBS
a) Wygładzanie: podział powierzchni NURBS na płatki wypukłe o odpowiednio małej krzywiźnie. Długi i skomplikowany proces – używany algorytm deCastlejau, dzielący płatki. [9]

b) Na podstawie tego podziału jest generowana hierarchia BV (Bounding Volumes), płatek po podziale nie jest przechowywany w pamięci. Cały podział był tylko po to aby uzyskać obszary na powierzchni NURBS (BV) w których mamy (ale nie na 100%) zagwarantowaną szybką zbieżność metody Newtona. W BV trzymamy współrzędne punktu startowego dla tego segmentu (u,v). W pamięci przechowujemy tylko oryginalny NURBS bez informacji o podziałach. Używane BV to AABB drzewa (Axis-aligned Bounding Boxes), OBB drzewa (Oriented Bounding Boxes) lub sfery (drogi koszt wyliczenia)

c) Promień przekształcamy na dwie powierzchnie, których płaszczyzną przecięcia jest promień. Będziemy szukać przecięć tych płaszczyzn i powierzchni NURBS. Najpierw będą wyznaczane przecięcia z BV a potem dopiero dla znalezionych punktów startowych metoda iteracyjna: np. Newtona.

d) Ostatecznie należy wyliczyć punkt metodą Newtona. W wielu przypadkach metoda jednak nie będzie zbieżna, należy ustawić maksymalną ilość kroków. Ogólnie całość działa bardzo wolno i pomimo użycia tylu narzędzi artefakty są najczęściej widoczne....



	1.	Metoda raytracingu powierzchni NURBS przez triangulację (stosowana w moim algorytmie)

W preprocesingu dzielę powierzchnię NURBS na siatkę d1 x d2 punktów. Liczby d1 i d2 są podane przez użytkownika w definicji powierzchni NURBS. Podział nie jest jednorodny ponieważ w wyniku triangulacji powierzchnia NURBS nie będzie gładka, ale będzie przybliżana płaskimi trójkątami. Dlatego na obszarach  położonych bliżej „krawędzi” powierzchni triangulacja jest „gęstsza” a na obszarach środkowych „rzadsza”.

W moim algorytmie jest możliwość wybrania metody triangulacji. Podając parametr –v liczba, możemy ustawić jak ma być dzielona powierzchnia NURBS. –v 1 ustawia na podział równomierny, większe wartości ustawiają podział gęstszy na brzegach, mniejsze gęstszy w środku.

	Współczynniki te można podać z wiersza poleceń programu, lub pominąć – wtedy zostaną użyte rozsądne wartości domyślne.

	W punktach podziału obliczana jest powierzchnia NURBS dla tych punktów. D1 dzieli u[0,1] na d1 punktów, podobnie d2 dzieliło v[0,1] na d2 punktów. Powstaje nam w ten sposób d1*d2 czworokątów. Każdy z tych czworokątów dzielę na dwa trójkąty. Ostatecznie mamy 2*d1*d2 trójkątów z powierzchni NURBS. Kolejne powierzchnie NURBS są dodawane do listy trójkątów na końcu. tj. jeżeli scena miała 1000 trójkątów i był na niej NURBS o 256 trójkątach to indeksy jego trójkątów będą z przedziału [1000,1256] (oryginalnej sceny: [0,999])


	2.4.1 Obliczanie normalnych powierzchni NURBS


	Aby tak ztriangulowana powierzchnia wyświetlała się dobrze, należy jeszcze obliczyć normalne we wszystkich wierzchołkach trójkątów. Normalna do powierzchni NURBS jest obliczana w przybliżeniu. Wektor normalny jest iloczynem wektorowym pochodnych w kierunku u i v, tj: n = Pv x Pu. Oto algorytm obliczający normalną do zadanej powierzchni NURBS:

Algorytm:

Normalna(u,v)

Wejście: u,v: punkt na powierzchni NURBS
Wyjście: x,y,z: współrzędne wektora normalnego

 EpsilonU = 0.000001
 EpsilonV = 0.000001

 jeżeli 1-2*EpsilonU <= u wtedy EpsilonU = -EpsilonU
 jeżeli 1-2*EpsilonV <= v wtedy EpsilonV = -EpsilonV

 (x1,y1,z1) = NURBS(u, v)
 (x2,y2,z2) = NURBS(u+EpsilonU, v)
 (x3,y3,z3) = NURBS(u, v+EpsilonV)

  dx1 = x2-x1
  dx2 = x3-x1
  dy1 = y2-y1
  dy2 = y3-y1
  dz1 = z2-z1
  dz2 = z3-z1

  x = dy1*dz2 - dz1*dy2
  y = dz1*dx2 - dx1*dz2
  z = dx1*dy2 - dy1*dx2

  długość = pierwiastek(x*x + y*y + z*z)
  (x,y,z) = (x,y,z) / długość

  jeżeli EpsilonU*EpsilonV > 0 to (x,y,z) = -(x,y,z)

	W ten sposób dla każdego trójkąta mamy obliczone 3 normalne w każdym z jego wierzchołków. Teraz aby obliczyć normalną w dowolnym punkcie trójkąta należy zastosować interpolację.


	2.4.2 Triangulacja, sposób podziału powierzchni


	W ten sposób podzieliłem powierzchnię NURBS na 2*d1*d2 trójkątów GPT (Generalized Phong Triangles). Używam trójkątów GPT aby uzyskać gładkie przejścia pomiędzy wierzchołkami i poszczególnymi trójkątami. Artefakty będą widoczne głównie na brzegu powierzchni NURBS (będzie wyglądała na łamaną), mogą się także ujawnić przy cieniach (ponieważ są rzucane na płaskie trójkąty). Przy dość gęstej triangulacji i odpowiednio rozmieszczonej artefakty albo słabo widoczne albo nie widać ich w ogóle. Zysk szybkości jest olbrzymi w porównaniu do obliczeń iteracyjnych metodą Newtona. Ztriangulowane powierzchnie NURBS są dodawane do sceny jako zwykłe listy trójkątów i tworzone jest dla nich drzewo lokalizacyjne AABB-Tree (Axis-Aligned Boundig Box Tree). Algorytmy tworzące drzewo AABB będą opisane poniżej.


	Wszystkie właściwości trójkąta GPT są interpolowane z jego wierzchołków. Podstawowe właściwości to: normalne, kolor, współrzędne tekstury, współczynniki: załamania, odbicia, pochłaniania, chropowatości i rozbłysku dla każdego koloru. Powierzchnia NURBS zawiera wszystkie te właściwości zapisane w każdym punkcie kontrolnym.

	2.4.3 Struktura trójkąta GPT

	Struktura opisująca trójkąt jest następująca:

Wierzchołek: x,y,z
Normalna: x,y,z
Materiał: kolor_r, kolor_g, kolor_b
	    odbicie_r, odbicie_g, odbicie_b
	    przezroczystość_r, przezroczystość_g, przezroczystość_b

Wartości odpowiadają kolejnym kolorom: czerwony, zielony, niebieski: r,g,b oraz właściwościom: rozproszenie (kolor), odbijalność (odbicie) i przezroczystość. Wszystkie wartości są z przedziału [0;1]. Powinny sumować się do 1, jeżeli się nie sumują to są skalowane np. dla koloru czerwonego (c-kolor, s-odbijalność(specular), t-przezroczystość(transparency): cr=1, tr=2, sr=1. Zostaną skalowane do: c=0.25, s=0.25, t=0.5. Poszczególne kolory r, g, b są niezależne.

Współrzędna tekstury: s,t
Współczynniki załamania mXX: pierwsza litera oznacza czy z ośrodka do powierzchni (D) lub z powierzchni do ośrodka (U). Druga odpowiada za kolor, mamy więc np. dla koloru czerwonego: Przechodząc z powierza do obszaru ograniczonego powierzchnią NURBS: używamy mDR, a wychodząc z obszaru ograniczonego powierzchnią: mUR. Ważne aby obszar ograniczony powierzchnią NURBS był zamknięty.

Współczynniki świetlistości.
	Odpowiadają za to jak „ostre” są rozbłyski światła na wypolerowanych powierzchniach osobno dla każdego koloru: czerwony, zielony i niebieski. Im większy ten współczynnik tym „ostrzejsze”. Odpowiada on za światło rozbłysku w modelu Phonga, jest to potęga do której podnosimy kosinus kąta między wektorem normalnym a wektorem skierowanym do obserwatora.

Trójkąt:
	Wierzchołek: a,b,c
	Normalna: na,nb,nc
	Materiał: ma, mb,mc
	Współrzędna Tekstury: ta,tb,tc
	Tekstura: t
	Powierzchnia: s
LiczbaRzeczywista: mUR, mUG, mUB, mDR, mDG, mDB
	LiczbaRzeczywista: wspólczynnik_świetlistościR
	LiczbaRzeczywista: wspólczynnik_świetlistościG
	LiczbaRzeczywista: wspólczynnik_świetlistościB
LiczbaCałkowita: indeks_NURBS (indeks pow. NURBS do której należy lub –1)


	2.4.4 Obliczanie powierzchni w której zawiera się trójkąt
Powierzchnia: A,B,C,D

	Definicja powierzchni może być dostarczona przez użytkownika, lub można podać 0,0,0,0 (nie reprezentujące żadnej powierzchni), wtedy program policzy powierzchnię sam. Opis powierzchni przez A,B,C,D: A,B,C to współczynniki normalnej tej powierzchni a D to funkcja odległości od środka układu współrzędnych (dla znormalizowanych A,B,C to jest odległość). Zapis nie jest jednoznaczny, można go skalować: tA,tB,tC,tD, dla t <> 0. powierzchnie są używane przy obliczaniu przecięć trójkątów i promieni.

	Algorytm obliczający powierzchnię dla zadanego trójkąta:

LiczPowierzchnie(Trojkat t)

Wejście: Trójkąt
Wyjście: Trójkąt z obliczoną powierzchnią w której się zawiera. (sA,sB,sC,sD)


 (x1,y1,z1)  =  B – A
 (x2,y2,z2)  =  C – A

 nx = y1*z2 - z1*y2
 ny = z1*x2 - x1*z2
 nz = x1*y2 - y1*x2

 długość = pierwiastek(nx*nx + ny*ny + nz*nz)
 (nx,ny,nz) = (nx,ny,nz) / długość

 sA = nx
 sB = ny
 sC = nz

 sD = -sA*A.x - sB*A.y - sC*A.z


	2.4.5 Obliczanie normalnych trójkąta na podstawie wierzchołków

	Wartości normalnych mogą zostać podane jako: (0,0,0), wtedy zostaną obliczone:

ObliczNormalne(t)

Wejście: Trójkąt: t
Wyjście: trójkąt z policzonymi normalnymi w wierzchołkach.


(x1,y1,z1) = B – A
(x2,y2,z2) = C – A

 nx = y1*z2 - z1*y2
 ny = z1*x2 - x1*z2
 nz = x1*y2 - y1*x2

 długość = pierwiastek(nx*nx + ny*ny + nz*nz);

(nx,ny,nz) = (nx,ny,nz) / długość

Na = Nb = Nc = (nx,ny,nz)

	2.4.6 Algorytmy interpolacji

	Wszystkie te wartości są interpolowane na powierzchni trójkąta. W szczególności interpolacja jest używana do obliczenia koordynatów tekstur na całym trójkącie. Są do wyboru dwa algorytmy interpolacji: bazowany na odległości od najbliższej krawędzi (wolniejszy i testy wykazały, że mniej odpowiedni) (1) oraz dwu-krokowy:  najpierw interpolacja na rzucie punktu przecięcia na jedną z krawędzi a potem interpolacja pomiędzy wynikiem pierwszego kroku a wierzchołkiem z którego wykonano rzut (szybszy i lepszy algorytm interpolacji dwuliniowej) (2). Oto algorytm pierwszy (1)

Algorytm pierwszy:


InterpolujAlg1(t,p)

Wejście: Trójkąt t, Wierzchołek p (należący do trójkąta)
Wyjście: wszystkie właściwości trójkąta interpolowane dla zadanego wierzchołka. Uwaga: funkcja zakłada, że wierzchołek należy do trójkąta, jest ona wywoływana gdy procedura Intersection (przecięcie) znajdzie punkt wspólny promienia i trójkąta.

 ab = Linia [B – A]
 bc = Linia [C – B]
 ca = Linia [C – A]
 ap = Linia [A – p]
 bp = Linia [B – p]
 cp = Linia [C – p]

 jeżeli pa = przecięcie_linii(bc, ap) nie istnieje to błąd
 jeżeli pb = przecięcie_linii(ca, bp) nie istnieje to błąd
 jeżeli pc = przecięcie_linii(ab, cp) nie istnieje to błąd


 la = odległość(pa, p)
 lb = odległość(pb, p)
 lc = odległość(pc, p)

 suma = la+lb+lc
 (la,lb,lc) = (la,lb,lc) / suma

 N = la*Na + lb*Nb + lc*Nc
 Unormuj(N)

(pozostałe właściwości podobnie)

Obsługa przypadków szczególnych (takich jak punkt p bardzo blisko któregoś wierzchołka – patrz rayslib.c get_normal_new).

Oto algorytm drugi (2). Domyślnie używany przez program, można jednak ustawić pierwszy.


Algorytm drugi

InterpolujAlg2(t,p)

Wejście: Trójkąt t, Wierzchołek p (należący do trójkąta)
Wyjście: Wszystkie właściwości trójkąta interpolowane dla zadanego wierzchołka. Uwaga: funkcja zakłada, że wierzchołek należy do trójkąta, jest ona wywoływana gdy procedura Intersection (przecięcie) znajdzie punkt wspólny promienia i trójkąta.


Normal tymczN

bc = Linia [C – B]
ap = Linia [P – A]

jeżeli d = przecięcie linii(bc, ap) nie istnieje to przypadek szczególny

 dlugoscBC = odległość(c, b)
 dlugoscBD = odległość(d, b)

 wsp = długośćBD / dlugośćBC

 tymczN = (1-wsp)*Nb + wsp*Nc
 unormuj(tymczN)

 długośćAD = odległość(d, a)
 długośćAP = odległość(p, a)

 wsp = lenAD / lenAP;

 N = (1-wsp)*Na + wsp*tymczN
 Unormuj(N)

	2.4.7 Algorytm przecięcia linii

Używając powyższych algorytmów interpolowane są wszystkie właściwości trójkątów GPT, efekt jest dobry: szczególnie dobrze do testów tego typu są interpolacje tekstury. Funkcja obliczająca przecięcia linii została starannie napisana aby uwzględniać wszystkie przypadki szczególne: równoległości, nieskończenie wiele punktów, linie pionowe, poziome, równoległe do z itd. Jej kod jest bardzo długi, przedstawię tutaj tylko kawałek pseudo-kodu obliczający typowy przypadek, reszta jest w rayslib.c: line_intersect. Metoda ta oblicza przecięcia promieni, struktura opisująca promień jest następująca:

Promień:
	Wierzchołek P jest to wierzchołek od którego biegnie promień
	Wektor d kierunek promienia
	LiczbaCałkowita r,c
		R – aktualna głębokość rekursji
		C – kolor: czerwony, zielony lub niebieski

PrzecięcieLinii(Promień p1, p2)

Wejście: dwa promienie
Wyjście: ilość punktów przecięcia: 0,1 lub nieskończenie wiele oraz ewentualny punkt przecięcia

(x1,y1,z1) = p1.P
(x2,y2,z2) = p2.P
(dx1,dy1,dz1) = p1.d
(dx2,dy2,dz2) = p2.d

(bez przypadków szczególnych!)

L = (dy1*(x2-x1)-dx1*(y2-y1)) / (dx1*dy2-dx2*dy1)
jeżeli dx1 <> 0 to K = (x2-x1+L*dx2)/dx1
wpp. jeżeli dy1 <> 0 to K = (y2-y1+L*dy2)/dy1
wpp K = 0
(x,y,z) = (x1+K*dx1, y1+K*dy1, z1+K*dz1)

(jedno rozwiązanie)

	2.5 Wady obliczania raytracingu przez triangulację

Wady przyjętego przeze mnie rozwiązania:

a) Widoczne artefakty (fragmenty łamanej) na krawędzi powierzchni NURBS przy zbyt „rzadkiej” triangulacji.

b) Źle rzucane cienie na powierzchnię (ponieważ są rzucane na płaskie trójkąty) – chyba największa wada. W algorytmie można wyłączyć/włączyć generowanie cieni.

c) Złe efekty gdy renderujemy powierzchnie pod dużymi kątami, szczególnie gdy powierzchnia jest do nas „bokiem”

d) Wszystkie te wady da się rozwiązać przez odpowiednio gęstą triangulację (proces jest zbieżny do rozwiązania dokładnego) ale wtedy narastają wymagania pamięciowe i czas preprocesingu (tworzenie drzewa AABB)

	2.6 Zalety przyjętego przeze mnie rozwiązania

a) Czas działania. Zysk kilku rzędów wielkości: 100x – 1000x.  Nie są konieczne kosztowne obliczenia numeryczne, a trójkąty można umieścić w strukturze lokalizacyjnej tak samo jak pozostałe trójkąty sceny. Szybkość i efekty są  zadowalające

b) Stabilność rozwiązania. Metoda Newtona boryka się z wieloma problemami zbieżności, przypadkami szczególnymi, wielokrotnymi rozwiązaniami, rozbieżnością itp. Jak już  wspomniałem powyżej wiele było prób napisania dokładnego algorytmu który byłby stabilny numerycznie i radził sobie ze wszystkimi przypadkami szczególnymi, większość tych prób zakończyła się niepowodzeniem. Moja metoda daje 100% stabilność i bardzo szybko generuje wynik, przy niewielkiej (relatywnie) stracie dokładności.

c) Zastosowanie trójkątów GPT i interpolacji pozwala usunąć prawie wszystkie problemy związane z tym, że trójkąty są płaskie

d) Ponieważ z powierzchni NURBS powstają listy trójkątów to możemy je przekształcać jak wszystkie inne obiekty w programie, np. przemieścić część, zmienić tekstury itd. Mamy także możliwość zapisania dobrego przybliżenia powierzchni NURBS jako trójkąty.




	3. Główne Algorytmy

	Przedstawię poniżej główne algorytmy używane przez mój program.

	3.1 Algorytmy przecięć

	Najważniejszym algorytmem Raytracingu jest niewątpliwie algorytm przecięcia: Intersection. Mój program przecina najpierw Hierarchię obiektów otaczających a dopiero na koniec przecina trójkąty. Tutaj opiszę algorytm najniższego poziomu tj. przecięcie promienia z trójkątem. W programie jest możliwość wyboru algorytmu przecinającego: stara wersja przecina po prostu promień ze wszystkimi trójkątami, a nowa używa BVH (Bounding Volume Hierarchy – AABB Tree). Poniżej przedstawiam algorytm wywoływany przez oba przecięcia na najniższym poziomie: przecięcie trójkąta i promienia:


PrzetnijTrojkąt(t,r)

Wejście: Trójkąt t, Promień r
Wyjście: Punkt przecięcia albo NULL

Powierzchnia s = t.S

wsp = s.A*r.d.x +  s.B*r.d.y + s.C*r.d.z
jeżeli wsp <> 0  to
t = -(s.A*r.P.x + s.B*r.P.y + s.C*r.P.z + s->D)/wsp;
p = r.P + r.d * t
jeżeli wierzcholek_w_trójkącie(t, p)
jeżeli po nieodpowiedniej stronie półprostej   (promienia) to zwróć NULL
wpp. Zwróć p
 zwróć NULL (promień równoległy do pow. Trójkąta)

	Algorytm testujący czy wierzchołek jest w trójkącie. Wykorzystuje sumę kątów między poszczególnymi bokami a punktem.

Wierzchołek_w_trójkącie(t,v)

Wejście: Trójkąt t, Punkt v
Wyjście: odpowiedź: tak lub nie

 Kąt = oblicz_kąt(w,a,b)
 Kąt = oblicz_kąt(w,c,c)
 Kąt = oblicz_kąt(w,c,a)
 Jeżeli kąt < 2*PI zwróć FALSE wpp. Zwróć TRUE


	3.2 Hierarchia  obiektów otaczających (BVH)

Krytycznym  czynnikiem decydującym o szybkości raytracingu nie jest procedura przecinająca (chociaż ma ona bardzo duże znaczenie), najważniejsza jest odpowiednia lokalizacja obliczeń. Załóżmy, że mamy scenę składającą się z 20000 trójkątów o rozdzielczości 1000 x 1000. promieni mamy więc 1000 x 1000 x 3 = 3M. Każdy z tych promieni musi przetestować przecięcie z każdym z trójkątów (gdy nie mamy BVH) tj. 3M x 20000 = 60G. Doliczyć należy jeszcze rekursję promieni. Niech każdy promień załamuje się i ugina, mamy wtedy 2 promienie z jednego na poziomie 1 rekursji. Na poziomie np. 5-tym mamy 2^5 = 32 promienie. Nasze 60G należy pomnożyć przez 32 mamy 1,9T. Zakładając jednak, że mamy strukturę lokalizacyjną która przechowuje w drzewie trójkąty mamy koszt przejścia przez drzewo z grubsza logarytmiczny (o podstawie 2). Więc Log2(20000)  < 15. Jednak nie zawsze możemy aż tak zredukować  koszt obliczeń, w tej sytuacji (wynika to z testów) możemy liczyć na 200-krotne zmniejszenie ilości obliczeń. Jest to zysk dosyć znaczący.

Drzewo AABB

Postanowiłem zaimplementować drzewo lokalizacji AABB-Tree. Nazwa z ang. Axis Aligned Bounding Box Tree – Drzewo wyrównanych do osi układu współrzędnych prostopadłościanów. Zastosowałem optymalizację drzewa użytą w bibliotece Opcode: mianowicie zlikwidowałem liście zawierające tylko jeden trójkąt, w najniższym poziomie drzewa (liściu) przechowywane są wskaźniki na dwa trójkąty, które ten liść otacza. Oznacza to konieczność wyszukania dwóch jak najbliższych trójkątów – co jest dosyć powszechne, gdyż modele są na ogół triangulowane przez siatkę prostokątną: każdy prostokąt to dwa trójkąty. Tak stworzone drzewo ma N-1 węzłów (dla N trójkątów), pełne drzewo miałoby 2N-1 węzłów.

Struktury danych BVH

Za przechowywanie danych BVH (hierarchii otaczającej) odpowiadają następujące struktury...

AABB: wyrównany do osi prostopadłościan:

Box:
	Punkt: min,max
	Wskaźnik: t1,t2 –trójkąty jeżeli jest to liść


Jeżeli na scenie jest nieparzysta ilość trójkątów to jeden box ma tylko jeden trójkąt.

Drzewo AABB:
Btree:
	BtreeWskaźnik: l,r
	Box: b
	Opiszę teraz algorytmy generowania, przecinania, wczytywania i zapisywania drzew AABB.

	3.2.1 Algorytm generacji drzewa AABB
	Są trzy możliwości wybrania algorytmu generującego drzewo: algorytm szybki, algorytm częściowy, algorytm pełny. Różnice polegają na sposobie minimalizacji wyboru boxa w każdym kroku tworzenia drzewa.
Algorytm preprocessingu sceny:

PreprocessingSceny(t[])

Wejście: tablica trójkątów
Wyjście: AABB-drzewo sceny

Dla i = 0 do ilość_trójkątów dodaj_na_liste(głowa, t[i])

Dopóki są trójkąty na liście „głowa”
	t1 = t2 = NULL
	minimalnaS = +nieskończoność
	s = minimalnaS
	hd = głowa
	minT = głowa
	Dopóki są trójkąty na liście „głowa”
		Jeżeli s = oblicz_powierzch_trójk(głowa) < minimalnaS
			MinimalnaS = s
			MinT = głowa
		Głowa.następny
      t1 = minT
      głowa = hd
      usuń_z_listy(głowa, minT)
      hd = głowa
      minT = głowa
      minS = +nieskończoność
      s = minS
	jeżeli są trójkąty na liście „głowa”
		Dopóki są trójkąty na liście “głowa”
			Jeżeli s = oblicz_box(głowa, t1) < minimalnaS
			MinimalnaS = s
				MinT = głowa
		Głowa.następny
      t2 = minT
      głowa = hd
      usuń_z_listy(głowa, minT)
      hd = głowa
dodaj_box(lista_boxów, t1, t2)
utwórz_drzewo_z_listy(lista_boxów)

	Omówię teraz kolejno wykorzystywane w algorytmie funkcje (poza tymi, które są oczywiste).

	Funkcja obliczająca powierzchnię trójkąta (używana aby zminimalizować wielkości powstających boxów)

Oblicz_powierzchnię_trójkąta(t)

Wejście: trójkąt
Wyjście: miara jego powierzchni (nie koniecznie dokładna powierzchnia, ale zachowuje relację <)

ca = C – A
ba = B - A

la = długość(ca)
lb = długość(ba)
 unormuj(ca)
 unormuj(ba)

 alfa = arcus_kosinus(iloczyn_skalarny(ca, ba))

 ZWRÓC 0.5 * sin(alfa) * la * lb

	Funkcja obliczająca objętość boxa, który będzie otaczał dwa trójkąty (liść po optymalizacji). Liść w drzewie zawiera dwa trójkąty. Funkcja jest obliczana przy wyszukiwaniu minimalnego boxa, gdy taki zostanie znaleziony to dwa trójkąty które otacza zostają dodane do listy boxów. Po przejściu listy trójkątów w ten sposób powstaje sufit(N/2) liści drzewa i należy dobudować drzewo w górę, robi to funkcja: utwórz_drzewo_z_listy.

oblicz_box(t1,t2)

Wejście: Trójkąty: t1 i t2
Wyjście: objętość boxa

 minx = sprawdź wszystkie x, wybierz minimalny
 maxx = sprawdź wszystkie x, wybierz maksymalny
 miny = sprawdź wszystkie y, wybierz minimalny
 maxy = sprawdź wszystkie y, wybierz maksymalny
 minz = sprawdź wszystkie z, wybierz minimalny
 maxz = sprawdź wszystkie z, wybierz maksymalny

 dx = maxx - minx
 dy = maxy - miny
 dz = maxz – minz
 ZWRÓC dx*dx + dy*dy + dz*dz

Funkcja dodaj_box działa analogicznie, poza tym, że obliczony box jest dodawany do listy.

Najważniejszą funkcją jest funkcja tworząca drzewo z wygenerowanej listy boxów-liści. Jest ona wywoływana na końcu algorytmu: utwórz_drzewo_z_listy(lista_boxów). To właśnie w tej funkcji są rozgraniczenia na pełną/częściową/zerową optymalizację drzewa. Podobnie jak poprzednio minimalizowaliśmy wielkość boxa otaczającego dwa trójkąty, teraz będziemy minimalizować wielkość boxa otaczającego dwa boxy. Funkcja obliczająca box otaczający dla dwóch boxów (a nie trójkątów która została już zaprezentowana) działa bardzo podobnie do opisanej już funkcji oblicz_box. Nazwa tej funkcji to oblicz_box2. Cały algorytm generuje rekurenycje drzewo boxów, każdy box otaczający jest minimalnym boxem otaczającym swoje dzieci. Korzeniem drzewa jest box otaczający całą scenę.

3.2.1.1 Oto algorytm generujący drzewo z listy:

Utwórz_drzewo_z_listy(lista)

Wejście: lista boxów-liści wygenerowana wcześniej
Wyjście: AABB-Tree zbudowane dla tej listy (btree)

 hd = lista
 głowa = NULL
 box1 = box2 = NULL
 ilość = 0
 Dopóki są elementy na „lista”
	tymcz.l = tymcz.r = NULL
    	 tymcz.b = lista
	dodaj_do_listy(głowa, tymcz)
 	lista.następny
 ilość = ilość_elementów(głowa)
 ogon = pobierz_ogon_listy(głowa)
 Dopóki ilość > 1
	Tutaj wybieramy jeden z algorytmów minimalizacji
	tymcz = znajdz_najblizsze_boxy(głowa, ogon, box1, box2)
      Usuń_z_listy(głowa, ogon, box1)
      Usuń_z_listy(głowa, ogon, box2)
	Dodaj_na_koniec_listy(głowa, ogon, tymcz)
	Ilość = ilość - 1
 btree = głowa

Ważną funkcją jest znajdz_najbliższe_boxy(), w zależności od tego jakiego typu funkcję wybierzemy na nią, otrzymamy algorytmy: bez minimalizacji, z minimalizacją częściową lub z pełną.

3.2.1.2 Algorytm minimalizacji pełnej:

Znajdz_najblizsze_boxy_algorytm_pełny(h, t, b1, b2)

Wejście: głowa i ogon listy boxów
Wyjście: wybrane boxy zwracane w b1 i b2, bt (pojedynczy węzeł drzewa)


 p1 = head
 minD = +nieskończoność
 b1 = b2 = = bl = br = NULL
 Dopóki są elementy na “p1”
	p2 = p1.następny
	Dopóki są elementy na “p2”
		Jeżeli d = oblicz_box2(p1, p2) < minD
			MinD = d
			bl = p1
			br = p2
   		p2.następny
    	p1 = p1.następny
bt.l = bl
bt.r = br
ZWRÓĆ bt

	3.2.1.3 Algorytm minimalizacji częściowej

Minimalizacja częściowa: zarówno p1 jak I p2 są ograniczone w przeszukiwaniu listy jakąś liczbą całkowitą, np 10% oznacza, że p1 i p2 przeszukają tylko 2000 z 20000 trójkątów. Można tak zrobić ponieważ, drzewo boxów jest posortowane po kolejnych rosnących boxach, w pierwszym etapie (gdy tworzono listę liści) w którym zawsze jest stosowana minimalizacja pełna. Wygenerowane w ten sposób drzewo nie jest optymalne i w czasie procesu raytracingu może generować więcej kolizji, ale z reguły nie jest to duży narzut czasu. Złożoność maleje z O(N^3) dla pełnej minimalizacji do O(N^2) dla braku minimalizacji. Minimalizacja częściowa jest gdzieś pomiędzy (zależy ile procent). Zalecane wartości procentowe (przełącznik –F) to 0-8 % - dają już bardzo dobre wyniki. Minimalizacja pełna jest zalecana wyłącznie dla małych scen (ilość trójkątów mniejsza niż 8000). Brak minimalizacji nie jest zalecany w żadnym przypadku, należy użyć dowolnie małej wartości np. 0.01%.

Algorytm częściowy:

Znajdz_najblizsze_boxy_alg_częściowy(h, t, b1, b2,k)

Wejście: głowa i ogon listy boxów, k mówi nam czy przerwać przeszukiwanie po i czy po j (k = 0 to po i, a k = 1 to po j)
Wyjście: wybrane boxy zwracane w b1 i b2, bt (pojedynczy węzeł drzewa)

 p1 = head
 minD = +nieskończoność
 b1 = b2 = = bl = br = NULL
 I = j = 0
 Dopóki są elementy na “p1”
	p2 = p1.następny
	Dopóki są elementy na “p2”
		Jeżeli d = oblicz_box2(p1, p2) < minD
			MinD = d
			bl = p1
			br = p2
		j = j + 1
		jeżeli k = 0 i j > ilosc_kroków to p2 = NULL
   		p2.następny
	i = i + 1
jeżeli k = 1 i i > ilosc_kroków to p1 = NULL
    	p1 = p1.następny
bt.l = bl
bt.r = br
ZWRÓĆ bt


3.2.1.4 Algorytm szybki

Minimalizacja zerowa (brak minimalizacji) przeszukuje tylko dla jednego trójkąta. Pierwszy trójkąt z listy jest ustalony, a dla drugiego przeszukujemy wszystkie. Jest to dosyć dobre rozwiązanie ponieważ w procesie tworzenia listy liści drzewa posortowaliśmy boxy. Metoda jest bardzo szybka (właściwie preprocesing dla tej metody jest prawie natychmiastowy) i daje dobre rezultaty, jednak czas obliczeń kolizji wzrasta. W algorytmie jest także możliwość przetasowania listy, drzewa nie będą wtedy tak bardzo niezrównoważone. Istnieją przypadki dla których algorytm szybki może utworzyć zdegenerowane (liniowe) drzewo.

Znajdz_najbliższe_boxy_algorytm_szybki(h, t, b1, b2)

Wejście: głowa i ogon listy boxów
Wyjście: wybrane boxy zwracane w b1 i b2, bt (pojedynczy węzeł drzewa)


 p1 = head
 minD = +nieskończoność
 b1 = b2 = = bl = br = NULL
 p2 = p1.następny
 Dopóki są elementy na “p2”
	Jeżeli d = oblicz_box2(p1, p2) < minD
		MinD = d
		bl = p1
		br = p2
   	p2.następny
bt.l = bl
bt.r = br
ZWRÓĆ bt

	Doświadczalnie okazało się, że nie należy przeprowadzać pełnej minimalizacji sceny o złożoności o(n^3), wystarczy przeprowadzić minimalizację częściową lub nawet nie używać minimalizacji o(n^2). Dla dużych scen czas preprocessingu spada z kilkunastu godzin do kilku minut, a czas generowania obrazu wzrasta nie więcej niż o rząd wielkości: najczęściej 1.1 do 4 krotnie. Istnieją jednak pewne sceny dla których drzewo wygenerowane algorytmem pomijającym minimalizację są bardzo nieoptymalne, w pesymistycznym przypadku wysokość drzewa może być nawet rzędu ilości trójkątów, dlatego należy nawet dla dużych scen zastosować minimalizację częściową, np 0.1%. Dla małych scen (<5000 trójkątów) warto stosować minimalizację pełną. Dla bardzo dużych scen składających się z wielu rozłącznych obiektów, należy wygenerować drzewa dla poszczególnych obiektów, a potem połączyć je programem BTREECONV. Doświadczalnie okazało się także, że poziomy rekursji powyżej 7-8 są nierozróżnialne dla człowieka (i na poziomie 8-bitowej rozdzielczości kanałów RGB), poza szczególnymi przypadkami np. kamera wewnątrz lustrzanego sześcianu z innym obiektem w centrum.


	3.2.1.5 Algorytm wokselowy

	Najszybszy dostępny algorytm wykorzystuje dwie dodatkowe techniki. Po pierwsze należy zauważyć, że najczęściej odwiedzanym miejscem w drzewie jest jego korzeń, im głębiej w drzewie tym szansa trafienia tam jest mniejsza. Dlatego stosowany jest algorytm “smart minimalize”. Jest to odmiana minimalizacji częściowej, im bliżej liści tym mniej trójkątów jest przeszukiwanych w celu znalezienia minimalnego boxa otaczającego, zysk szybkości jest duży ponieważ liści jest w drzewie najwięcej a tam szukamy najmniej, a im wyżej w górę drzewa tym lepszą minimalizację stosujemy, ale takich elementów jest dużo mniej niż liści. Zysk na czasie jest bardzo duży (rzędu kilkukrotny), strata jakości drzewa niewielka (mierzona w ilości przecięć drzewa przez promienie, oraz ilości przecięć trójkątów w stosunku do przecięć boxów). Strata jakości drzewa wynosi od kilku do maksymalnie kilkudziesięciu procent (nie więcej niż 30%-40%). Aby wybrać ten algorytm należy zastosować opcję -F 200.

	Drugą techniką jest zastosowanie wokseli. Algorytm posiada dwa parametry, K-ilość podziałów oraz N-przy jakiej ilości trójkątów stosować podział. Domyślne wartości to K=2 i N=10000. Ustawienie N=10000 oznacza ze woksel będzie dzielony na K*K*K wokseli gdy ilość trójkątów w nim przekroczy 7500. Najlepsze wyniki uzyskuje się ustawiając K na 2 (oct-tree przestrzenie) Biorąc np. Scenę o 100000 trójkątów, algorytm działa rekurencyjne: Dzieli box otaczający całą scenę na 2*2*2 = 8 wokseli (o równych objętościach). Dla każdego z wokseli wykonuje procedurę rekurencyjne, aż do momentu gdy w wokselu będzie mniej niż 10000 trójkątów, gdy tak się stanie zastosuje dla nich algorytm “smart minimalize” tworząc box wynikowy. W fazie łączenia wszystkie boxy zostaną potraktowane jak dane wejściowe do algorytmu minimalizacji pełnej lub “smart” gdy będzie ich bardzo dużo. Dla bardzo dużych N (większych od ilości trójkątów na scenie) algorytm. Dla dużych K algorytm tworzy bardzo nieoptymalne drzewa, ponieważ niektóre woksele mają bardzo mało trójkątów juz na pierwszym poziomie rekursji, a inne są wielokrotnie. A potem w fazie łączenia są traktowane jednakowo, więc ścieżka do niektórych trójkątów jest bardzo długa, a do innych bardzo krótka. Polecam stosowanie K=2,3...6 (chociaż program pozwala stosować nawet K=32). Co do N to im większe tym lepiej, ale przedłuża to czas preprocesingu. Dla małych N algorytm ma złożoność liniową O(N*N*ilosc_trójkątów), polecane jest takie N dla którego pełna minimalizacja nie jest jeszcze bardzo długa, np N=5000. Zysk szybkości jest następujący: Weźmy scenę z 100000 trójkątów. Pełna minimalizacja to O(N^3) --> 10^15, Minimalizacja szybka to O(N^2) = 10^10, minimalizacja przez woksele nie większe niż 5000 to: Ponieważ 100K / 5K = 20 to w najlepszym przypadku mamy 20 wokseli, wszystkie po 5000 trójkątów, ale tak się nie zdaża, załóżmy że otrzymamy około 30-50 wokseli. Obliczenie jednego maksymalnie dużego woksela to 5000^3 (pełna minimalizacja) lub 5000^2 (najszybsza minimalizacja) czyli od 2.5*10^7 do 1.25 * 10^11. należy jeszcze zminimalizować woksele wynikowe (jest ich mniej niż 100, pełna minimalizacja to 100^3 = 10^6 co nie ma znaczenia przy koszcie ich obliczenia od 2.5*10^7 do 1.25*10^11). Niech ilość wokseli to 50. Ogólny koszt to od około 1.25*10^9 do 6.125*10^12. porównajmy to z kosztami algorytmu bez wokseli:
	Woksele (pełna minimalizacja) / Pełna Minimalizacja  ~= 150
	Woksele (smart minimalize)    / Pełna minimalizacja  ~=10000
	Woksele (fast minimalize)       / Pełna minimalizacja ~= 1,000,000
	Zalecanym i najlepiej spisującym się algorytmem jest smart minimalize używającym wokseli. Parametry algorytmu wybiera się poprzez -H “K N”, np -H “2 5000”, oznacza wybranie podziałow na 2*2*2=8 wokseli w każdym kroku rekurencyjnm, z warunkiem stopu rekurencji: ilość trójkatów w wokselu < 10000. Aby wymusić dodatkowo algorytm “full minimalize” należy dodać -F -100, dla “smart” -F -200 (ale to jest domyslna opcja więc nie trzeba jej podawać)

Algorytm:

Polącz_drzewa_wokselowe(l, n)
Wejscie: lista boxów: l, ilość boźów: n
Wyjście: połączone boxy w drzewo

Dopuki n > 1
 	b = znajdz_najblizsze_boxy(l, b1, b2)
	unun_z_listy(l ,b1)
	unun_z_listy(l ,b2)
	dodaj_do_listy(l, b)
     n = n -1


Oblicz_woksel(w, N, K, tab)
Wejście: Woksel: w, liczby K,N (opisane w 2.1.4.5)
		tab tablica znaczników czy trójkąt jest używany
Wyjście: przetworzony woksel


Jeżeli woksel.ilość_trójkątów > N
	n = N^3
	dla i=0 do N
		dx = maxx – minx
		dx = maxy – miny
		dx = maxz – minz
		iz = i % K
		iy = (i / K) % K
		ix = i / (K^2)
		w.woksel[i].minx = minx + ix/K * dx
		w.woksel[i].miny = miny + iy/K * dy
		w.woksel[i].minz = minz + iz/K * dz
		w.woksel[i].maxx = maxx + (ix+1)/K * dx
		w.woksel[i].maxy = maxy + (iy+1)/K * dy
		w.woksel[i].maxz = maxz + (iz+1)/K * dz

	     w.woksel[i].ilość_trójkatow = 0
	dopuki są trójkąty na liście l
		jezeli tab[l.t.idx] == 1 to pomiń trójkąt (już dodany)
		tab[l.t.idx] = 1

		x = (l.t.a.x + t.t.b.x + l.t.c.x) / 3
		y = (l.t.a.y + t.t.b.y + l.t.c.y) / 3
		z = (l.t.a.z + t.t.b.z + l.t.c.z) / 3

       x = (t->a.x + t->b.x + t->c.x) / 3.;
       y = (t->a.x + t->b.x + t->c.x) / 3.;
       z = (t->a.x + t->b.x + t->c.x) / 3.;

	  ix = (x-minx) / (maxx–minx) * K
	  iy = (y-miny) / (maxy–miny) * K
	  iz = (z-minz) / (maxz–minz) * K
       idx = K^2 * ix + K * iy + iz;
	  dodaj_do_listy(w.woksel[idx].lista_trójkątów, l.t)
       w.woksel[idx].ilość_trójkątów ++;
       pobierz następny trójkąt
    dla i=0 do N
		jezeli w.wokse[i].ilosc_trójkątów > 0  to
			ObliczWoksel(w.waoksel[i], N, K, tab)
wpp. Ilość trójkątów w wokselu < N
	standardowe obliczenie drzewa którymś z algorytmów 	minimalizacji
	dodaj box do listy

Oblicz_scene_uzywając_wokseli(t)
Wejście: lista trójkątów: t
Wyjście: obliczona scena

Woksel w		woksel glowny (obejmuje całą scenę)
w.box.minx = 1e10
w.box.miny = 1e10
w.box.minz = 1e10
w.box.maxx = -1e10
w.box.maxy = -1e10
w.box.maxz = -1e10

dla i=0 do ilość_trójkątów tablica[i] = 0	czy już używany?
dla i=0 do ilość_trójkątów
	dodaj_do_listy(w.lista_trójkątów, t[i])
	wyszukuj maksymalne i minimalne x,y,z
	ze wszystkich trójkątów
Przypisz te maksymalne i minimalne wartości do w.box
w.ilość_trójkątów = ilość_trójkątów
 v_root.ntbox = nTriangles;
 lista_wokseli.pusta()
 ObliczWoksel(w, N, K, tablica)
 Polącz_Drzewa_Wokselowe(lista_wokseli, ilość_trójkątów)


	3.2.1.6 Sąsiedztwo trójkątów

	Ponieważ w jeden piksel są wystrzelana 3 promienie (czerwony, zielony i niebieski), przecięcia ze sceną są obliczane tylko raz i zapisywane w drzewie indeksowym (struktura drzewa jest taka jak struktura wystrzelonego promienia, tj. Ma gałęzie odpowiadające za odbicie, załamanie i obliczenie cienia, które są używane wtedy gdy promień tworzy w trakcie obliczeń wyżej wymienione promienie: odbity, załamany lub cień) Algorytm posiada parametr M, określający co ile pikseli obliczać przecięcie. Gdy obliczamy piksel o numerze równym wielokrotności M, to do obliczeń używamy algorytmu przecinającego całe AABB-drzewo, wpp ustawiamy indeksy trójkątów sprawdzanych przez drzewo indeksowe na -1 (oznaczając w ten sposób, że indeksy te należy policzyć od nowa). Obliczając przecięcie, sprawdzamy rekurencyjne indeksy trójkątów w drzewie, jeżeli indeks jest dodatni to testujemy przecięcie z trójkątem o tym indeksie, jeżeli przecięcie jest to oznaczamy je jako aktualne przecięcie i obliczamy następny piksel (uwaga – jest to uproszczenie, inne przecięcie też może istnieć i być bliżej źródła promienia, niż to obliczone, z tego powodu M powinno być równe 1, ponieważ obliczenia pozostałych kolorów dla tego samego piksela nie spowodują powstania przecięć bliżej obserwatora, ale w kolejnym pikselu taka sytuacja może juz zaistnieć, dlatego domyślną wartością M jest 1). należy także zwrócić uwagę na załamanie światła dla M=1, ponieważ załamane promienie mogą mieć różne kierunki, to obliczenia przecięcia dla kolejnych kolorów w tym samym pikselu, także mogą dać niepoprawne wyniki, jednak przekłamania są niewielkie i zostały przeze mnie zaakceptowane, aby wyłączyć obliczanie z użyciem sąsiedztwa, należy ustawić M=0). Jeżeli przy testowaniu przecięcia z indeksem w drzewie da wynik negatywny, to wartość tego indeksu jest ustawiana na trójkąt z którym jest przecięcie, obliczany za pomocą algorytmu standardowego używającego drzewa AABB. Na tym etapie można także wybrać algorytm wykrywający krawędzie trójkątów (wtedy gdy przechodzimy do kolejnego trójkąta następuje zmiana indeksu w drzewie, wystarczy zwrócić wtedy wynik testu jako TRUE, i w krawędziach otrzymamy kolor tła. Ustawianie dużych M (np > 10) daje algorytm znacznie szybszy, ale jego dokładność jest wysoce niezadowalająca!


	3.2.2 Zapis i odczyt drzewa AABB

Po wykonaniu wszystkich opisanych powyżej kroków mamy wygenerowane AABB drzewo. Tak wygenerowane drzewo możemy następnie zapisać w pliku, dodając do każdego węzła jego numer (np. przechodząc drzewo pre-order lub dowolnie) i zamieniając wszystkie wskaźniki na indeksy odpowiednich węzłów. Należy także zapisać indeks korzenia drzewa. Podobnie możemy wczytać takie drzewo, generując tablicę węzłów i czytając indeksy zamieniać je na wskaźniki na odpowiednie elementy tablicy. W pliku mamy zapisany indeks korzenia. Algorytmy zapisu i odczytu wyglądają następująco:


3.2.2.1 Zapis drzewa AABB

ZapiszDrzewo(t)

Wejście: AABB-drzewo t
Wyjście: plik z zapisanym drzewem

 policz_indeksy()

 n = idx
 utwórz_tablice_powiązań(tab_idx)
 strcpy(fn, scenef);
 zapisz n
 idx = 0
 dla i = 0 do n jeżeli tab_idx[i] = korzeń zapisz i
 Zapisz_drzewo_do_pliku(korzeń, tab_idx,n)
 zamknij plik

 policz_indeksy()

 Btree tymcz
 idx = 0
 tymcz = korzeń
 przejdź_drzewo(tymcz)
 tutaj mamy policzony indeks w idx


 przejdź_drzewo(Btree t)

 jeżeli jest prawa gałąź to przejdź_drzewo(t.r)
 jeżeli jest lewa gałąź to przejdź_drzewo(t.l)
 idx = idx + 1


 utwórz_tablice_powiązań(tablica tab)

 Btree tymcz
 Idx = 0
 Tymcz = korzeń
 Przejdź_drzewo_uzupełniając_tablice(tymcz, tab)

 Przejdź_drzewo_uzupełniając_tablice(Btree t, Tablica tab)

 jeżeli jest prawa gałąź to przejdź_drzewo_uzupełniając_tablice(t.r, tab)
 jeżeli jest lewa gałąź to przejdź_drzewo_uzupełniając_tablice (t.l, tab)
 tab[idx] = t
 idx = idx + 1


Zapisz_drzewo_do_pliku(t, tab, n)

Wejście: drzewo (t), tablica indeksów (tab), wielkość tablicy (n)
Wyjście: plik z zapisanym drzewem AABB

jeżeli jest prawa gałąź to Zapisz_drzewo_do_pliku(t.r, tab, n)
jeżeli jest lewa gałąź to Zapisz_drzewo_do_pliku(t.l, tab, n)

zapisz indeks aktualnego węzła: idx
jeżeli jest prawa gałąź
	dla i = 0 do n jeżeli tab[i] = t->r to zapisz indeks prawy: i
	jeżeli nie znaleziono to zapisz –1
jeżeli jest lewa gałąź
	dla i = 0 do n jeżeli tab[i] = t->l to zapisz indeks lewy: i
	jeżeli nie znaleziono to zapisz -1
zapisz box w pliku (to jest oczywiste)
idx = idx + 1
	3.2.2.2 Odczyt drzewa AABB

Algorytm odczytujący zapisane AABB drzewo w pliku:

OdczytajDrzewo(t)

Wejście: plik z zapisanym drzewem AABB
Wyjście: drzewo AABB

odczytaj n (ilość elementów)
odczytaj root (indeks korzenia)
Jeżeli root <> n-1 to jest błąd

Utwórz tablicę n węzłów btree
Od i = 0 do n
	Odczytaj indeks i-tego elementu
	Odczytaj indeks prawej i lewej gałęzi
	Jeżeli prawy lub lewy indeks >= 0 to przypisz
		Btree[i].r = btree[indeks r]
		Btree[i].l = btree[indeks l]
	Odczytaj box
 btree = btree[root]
 Zamknij plik.

	3.2.3 Algorytm przecięcia (Intersection)

	Potrafimy już wygenerować oraz zapisać/odczytać drzewo. Opiszę teraz najważniejszy algorytm dla którego AABB drzewo zostało stworzone: mianowicie obliczanie przecięć promieni ze sceną. Idea algorytmu jest następująca:


	a) przecinamy promień z korzeniem drzewa (głównym największym boxem). Jeżeli nie ma przecięcia to zakańczamy algorytm, wpp. testujemy tą samą funkcją przecięcie z lewą i prawą gałęzią

b) jeżeli box jest liściem i wystąpiło przecięcie z nim, to przechodzimy do testowania przecięć z trójkątami.

	c) cały algorytm jest rekurencyjny i eliminuje niepotrzebne obliczenia przecięć z wieloma obiektami sceny, zysk często bywa 100 krotny lub 	nawet więcej, w zależności od wielkości sceny. Ilość testowanych przecięć jest rzędu O(Log2(N))

	Algorytm przecięcia może być wybrany spośród dwóch rodzajów: wykorzystujący AABB drzewo oraz nie (ten był napisany wcześniej na początku implementacji). Pomimo, że możliwość wyboru algorytmu istnieje, drugi algorytm jest stanowczo nie zalecany ze względu na jego bardzo powolne działanie




Algorytm intersekcji (przecięcia):

Intersection(t, r, v ,idx)


Wejście: Trójkąt: t
	    Promień: r
	    Punkt: v (ewentualny punkt przecięcia)
	    LiczbaCałkowita: idx (indeks kolidującego trójkąta)
Wyjście: Punkt: v, LiczbaCałkowita: idx

 rozwiązania = 0

 przecięcie_drzewa(głowa, r, btree)

 sol = rozwiązania;
 jeżeli sol = 1 to
	v = głowa.P
	idx = głowa.idx
	KONIEC
 Jeżeli sol < 1 to KONIEC

 min_odległość = +nieskończoność
 wielozn = NULL
 hd = głowa
 Dopóki są elementy na liście „głowa”
	Jeżeli odl = odleglosc(r.P, głowa.P) < min_odległość
       	Min_odległość = odl
		Idx = głowa.idx
		V = głowa.P
	Głowa.następny
 głowa = hd
 KONIEC


Przecięcie_drzewa(lista, r, tree)

Wejście: Lista lista (na niej zwracamy wszystkie przecięcia promienia z trójkątami), na początku ostatnie_rozwiązanie jest puste i lista jest pusta
	   Promień: r (promień przecinający)
	   Btree: tree (drzewo z którym przecinamy)
	   ostatnie rozwiązanie: Wierzchołek, na początku pusty, zmienna globalna
Wyjście: Lista punktów przecięcia z indeksami do odpowiednich trójkątów

Jeżeli nie ma gałęzi (ani prawe ani lewej)
	Jeżeli ret = PrzetnijTrójkąt(tree.b.t1, r)
		Lista.dodaj(ret, tree.b.t1.idx)
		Rozwiązania = rozwiązania + 1
		jeżeli bliżej początku promienia to ostatnie rozwiązanie = ret
	Jeżeli jest drugi trójkąt i ret = PrzetnijTrójkąt(tree.b.t2, r)
		Lista.dodaj(ret, tree.b.t2.idx)
		Rozwiązania = rozwiązania + 1
		jeżeli bliżej początku promienia to ostatnie rozwiązanie = ret
b = tree.r.box (sprawdzanie prawej gałęzi)
Przypadek_prawy_min_z:
K = (b.minz – r.P.z) / r.d.z (sprawdzenie po przedniej ścianie Z)
Jeżeli K > 0
	jeżeli są już przecięcia to:
		jeżeli (r.d.z > 0) L = (b.minz -ostatnie_przecięcie.z) / r.d.z
		jeżeli (r.d.z < 0) L = (b.maxz -ostatnie_przecięcie.z) / r.d.z
		jeżeli L > 0 to przeskocz do przypadek_prawy_min_y
x = r.P.x + K * r.d.x;
y = r.P.y + K * r.d.y;
jeżeli x >= b.minx i x <= b.maxx i y >= b.miny i y <= b.maxy
		Przecięcie_drzewa(lista, r, tree.r)
		Rozpatrz lewą gałąź
wpp. Jeżeli K < 0 to przeskocz do przypadek_prawy_min_y (pomiń max_z)
Przypadek_prawy_max_z:
K = (b.maxz – r.P.z) / r.d.z (sprawdzenie po tylnej ścianie Z)
Jeżeli K > 0
	jeżeli są już przecięcia to:
		jeżeli (r.d.z > 0) L = (b.minz -ostatnie_przecięcie.z) / r.d.z
		jeżeli (r.d.z < 0) L = (b.maxz -ostatnie_przecięcie.z) / r.d.z
		jeżeli L > 0 to przeskocz do przypadek_prawy_min_y
x = r.P.x + K * r.d.x;
y = r.P.y + K * r.d.y;
jeżeli x >= b.minx i x <= b.maxx i y >= b.miny i y <= b.maxy
		Przecięcie_drzewa(lista, r, tree.r)
		Rozpatrz lewą gałąź
Przypadek_prawy_min_y:
Podobnie z Y i X a potem lewa gałąź też Z,Y,X
		(...)
Prawdziwa implementacja testuje jeszcze 4 ściany (bo sprawdziliśmy tylko tylną i przednią ścianę) oraz także lewą gałąź. Tutaj pokazany został tylko jeden przypadek. Funkcja przecinania z trójkątem została pokazana wcześniej.


	3.3 Algorytm anty-aliasingu i rzucania cieni

	Zaimplementowany algorytm anty-aliasingu jest bardzo prosty: obliczenia prowadzone są w 2X większej rozdzielczości a następnie uśredniane są 4 piksele (macierz 2x2). Powoduje to zlikwidowanie niektórych nieprzyjemnych efektów: ostre krawędzie, aliasing tekstur itp.
Ponieważ promienie czerwony, zielony i niebieski są obliczane oddzielnie, więc jest możliwość uzyskania efektów rozszczepienia światła.
Aby zlikwidować efekt głębokich, czarnych cieni zastosowałem metodę półcienia. Można określić w programie minimalny i maksymalny poziom cienia, elementy półprzezroczyste rzucają półcienie, w zależności od ich współczynnika przepuszczania światła. Współczynnik ten może być różny dla różnych barw, więc możliwe jest otrzymanie kolorowych cieni. Problem z cieniami jest taki, że brane są pod uwagę tylko pierwsze obiekty na drodze: obiekt-światło. Jeżeli obiekt ten będzie prawie przezroczysty, to rzucony przez niego cień będzie bardzo mało intensywny. Jeżeli po drodze do światła dalej będzie zupełnie nieprzezroczysty obiekt to cień rzucany powinien być zupełny, ale mój algorytm obliczy ten słaby cień od pierwszego napotkanego obiektu. Jak widać mój algorytm nie jest doskonały ale dla dokładnego wyznaczania cieni lepiej używać innego algorytmu niż raytracing: np. „forward raytracing” opisany we wstępie 1. opcję generowania cieni można włączać/wyłączać. Można także ustawiać wartości minimalnego i maksymalnego zacienienia, jeżeli maksymalne zacienienie jest ustawione np. na 0.7, to współczynnik rzucania cienia przez obiekty jest mnożony przez 0.7.


3.4 Efekt Fresnela

Zaimplementowany został także efekt Fresnela. Zgodnie z prawem Fresnela gdy promień światła załamuje się pod dużym kątem w stosunku do normalnej materiału to „duża część” tego promienia zostaje odbita. Im większy jest kąt z normalną tym więcej światła się odbije a mniej przeniknie. Zaś gdy kąt padania jest bardzo zbliżony do normalnej to przezroczystość dąży do maksimum
W moim algorytmie założyłem, że dla promienia równoległego do normalnej mamy określone warunki w trójkącie (tj. współczynniki odnoszą się do kierunku normalnego). A wraz ze wzrostem kąta pomiędzy normalną a promieniem współczynnik załamania maleje na korzyść współczynnika odbicia (suma jest stała). Można zdefiniować z jaką potęgom jest zachowana ta proporcjonalność za pomocą współczynnika –k. Podanie –k 0 wyłącza efekt Fresnela, a np. –k 2 oznacza proporcjonalność kwadratową ze względu na iloczyn skalarny unormowanych wektorów. Np. równy 1 da w wyniku zawsze 1, ale np. 0.5 da: 1 dla –k 0 ale da 0.25 dla –k 2. Domyślna wartość to 0.625


	3.5 Algorytm wyznaczania koloru tła

	Zakładam że jeżeli promień nie natrafi na żaden obiekt, to jego kolor jest obliczany na podstawie tła. Tło może być zadane konkretnym kolorem przez użytkownika, lub można podać teksturę jako tło. Kolor tła jest podawany przez użytkownika i może być albo jednolity, albo obliczany na podstawie kierunku (8 kolorów na podstawie kierunków x,y,z: +,-). W przypadku tekstury współrzędne są obliczane na podstawie kierunku promienia zgodnie z następującym algorytmem:

	Algorytm obliczania koloru tła:

Kolor_tła(r)

Wejście: Promień r
Wyjście: wartości r,g,b koloru
Jeżeli nie zdefiniowano tekstury podłoża
	Jeżeli tryb jednokolorowy to
	r = TłoR
	g= TłoG
	b = TłoB
Jeżeli tryb 8 kolorów to
	jeżeli r.d.x > 0 to r = TłoR wpp r = 1 - TłoR
	jeżeli r.d.y > 0 to g = TłoG wpp g = 1 – TłoG
	jeżeli r.d.z > 0 to b = TłoB wpp b = 1 – TłoB
wpp.
	v = r.d
		unormuj(r.d)
	x = |v.x|
		y = |v.y|
	z = |v.z
|
	r = r + TłoR
	g =g + TłoG
	b = b + TłoB

	jeżeli r > 1 to r = 2 – r
	jeżeli g > 1 to g = 2 – g
	jeżeli b > 1 to b = 2 - b
wpp.
	v = r.d
	unormuj(r.d)
	x = |v.x|
	y = |v.y|
	z = |v.z|
	jeżeli x >=y i x >= z
	       xc = y
      	 yc = z
      inaczej jeśli y >=x i y >= z
		 xc = x
		 yc = z
	wpp
		 xc = x
		 yc = y
	pobierz_teksture(xc,yc,red,green,blue)
	zwróć (red,green,blue)


Nie jest to jedyny algorytm wyboru tekstury na podstawie kierunku, ale daje dość dobre efekty i dlatego go zastosowałem. Zalecenie: Tekstura powinna mieć rozmiary co najmniej 2x większe niż ekran i być samopowtarzalna.

	3.6 Algorytm raytracingu

Przedstawię teraz pseudo-kod algorytmu raytracingu:

Raytrace(s)

Wejście: Screen s (ekran)
Wyjście: obraz wynikowy

 odczytaj_scenę(trójkąty)

 btree = NULL
 jeżeli nie ma preprocesingu zapisanego PreprocesingSceny(trójkąty)
 wpp. OdczytajDrzewo(trójkąty)

 ZapiszDrzewo(trójkąty)

 r.P = observer
 r.r = 0

 dla i = 0 do s.x
    	r.d.x = i – s.x/2
	dla j = 0 do s.y
    		r.d.y = j – s.y/2
      oblicz_kolor(s, r, trójkąty, i, j)

 zapisz_obraz(s)


oblicz_kolor(s, r, t, x, y)

Wejście:  Screen: s
		Ray:	r
		Trójkąty[] t
		LiczbaCałkowita x,y
Wyjście: kolor piksela (x,y)

Jezeli nie ma intersection(t, r, v, idx)
       Kolor_tła(r, red, green, blue)
       Ustaw_kolor(s, x, y, red, green, blue)
wpp
	  InterpolujALg2(t, v)
       r.c = RED
       c_rec = 0
       red = kolor_rekurencyjny(r, t, v, n, idx, 1)
	  (tak samo green, blue)
       Ustaw_kolor(s, x, y, red, green, blue)


	3.6.1 Główna funkcja raytracingu – kolor_rekurencyjny:

	Maksymalna głębokość rekursji może być określona przez użytkownika, lub zostanie użyta domyślna: 6. Promień nie jest także obliczany gdy jego wkład w kolor jest mniejszy niż 1/256, ponieważ mamy tylko po 8 bitów na kolor na wyświetlaczu i taki piksel i tak nic by mnie mógł wnieść do wynikowego koloru.


Kolor_rekurencyjny(r, trójkąty, v, n, idx, f)

Wejście: Ray r
	   Triangle[] trójkąty
	   Punkt v (dla którego obliczamy kolor)
	   Normalna n (normalna w tym punkcie)
   LiczbaRzeczywista f (mnożnik – jak istotny jest ten promień w tym wywołaniu rekurencyjnym)
Wyjście: Zwraca wartość danego koloru (r.c) jako liczba rzeczywista


 t = trójkąty[id]

 Jeżeli r.r > maksymalna_rekursja
	Kolor_tła(r, red, green, blue)
	jeżeli r.c = RED ZWRÓĆ red
	(green,blue...)

 rv = 0
 Jeżeli źródło światła jest punktowe to tr = light – v
Wpp (wektorowe) tr = -light

ObsV = r.P - v
Jeżeli długość(obsV) < 1e-6 ZWRÓC 0
Jeżeli długość(tr) < 1e-6 ZWRÓC 0


 Ray li

 li.P = v
 li.d = tr
 li.r = r.r + 1

 unormuj(tr)
 unormuj(obsV)

 shadow = 1

	--cienie

 jeżeli jest intersection(trójkąty, li, nv, idx) i włączone obliczanie cieni
	jeżeli źródło światła jest punktowe
		l1 = odległość(v, light)
	wpp l1 = +nieskończoność
	l2 = odległość(v, nv)
	jeżeli l2 < l1
      	InterpolujAlg2(trójkąty[idx], nv)
		Jeżeli r.c = RED to shadow = (1-maxshadow) * nR.transparency
		(green,blue...)
		Jeżeli shadow > minshadow to shadow = minshadow

 Jeżeli shadow < maxshadow to shadow = maxshadow

	--światło rozproszone

 ca = iloczyn_skalarny(tr, n)

Jęeli r.c = RED to ca = ca * light.red_color
(g,b,…)

 Jeżeli r.c = RED to fact = R.c
 (g,b...)

	--tekstury

 Jeżeli trójkąt ma teksturę
	Pobierz_teksturę(t.tex, tC, red, green, blue)
	Jeżeli r.c = RED to fact = fact * red
	(g,b..)

 a2 = 2 * iloczyn_skalarny(obsV, n)

	--efekt Fresnela

 fresnel_t = |a2/2.|^fresnel_power
 fresnel_s = 1. - fresnel_t


 vi.P = v
 vi.d.x = a2*n - obsV
 vi.r = r.r + 1


 unormuj(vi.d)

 Jeżeli ca < ambient to ca = ambient

 rv = rv + fact * ca * shadow

	--rozbłysk Phonga
jeżeli r.c = RED to cca = t.caR
(...)

 jeżeli cca > 0
ca = potęga(iloczyn_skalarny(vi.d, tr), cca)
Jęeli r.c = RED to ca = ca * light.red_color
(g,b,…)
jeżeli ca > 0 to rv = rv + shadow * ca^ t.ca

 jeżeli rv > 1 to rv = 1
 Jeżeli r.c = RED to fact = fact = R.s + fresnel_s * R.t
 (g,b..)

 nf = f * fact

	--odbicie
Jeżeli nf > minimum_ray i rv < 1-minimum_ray to

    Jeżeli nie ma intersection(trójkąty, vi, nv, idx)

        	Kolor_tła(vi, &b_r, red, green, blue)
        	Jeżeli r.c = RED to fact2 = red
	  	(g,b,...)
		rv = rv + nf * fact2
	wpp.
		InterpolujAlg2(trójkąty[idx], nv)
		add = nf*rekurencyjny_kolor(vi, trojkąty, nv, nn, idx, nf)
	 	rv = rv + add
 Jeżeli r.c = RED
	fact = fresnel_t * R.t
	mU = t.mUR
	mD = t.mDR
 (g,b,...)

 nf = f * fact

		--załamanie

Jeżeli nf > minimum_ray i rv < 1-minimum_ray to
	ti.P = v
      a2 = iloczyn_skalarny(obsV, n)
	jeżeli a2 = 0 to skocz do ret_val
	jeżeli a2 > 0
    		--załamanie z ośrodka mU do mD
	       arg = ((sinus(arkus_kosinus(a2)))*mU)/mD
	       jeżeli arg < 0.9999 to beta = arcus_sinus(arg)
      	 wpp beta = 1.5705
        	 arg = kosinus(beta)

		 pn = -n * arg

      	 arg = mU/mD
		 pv = (a2*n – obsV)*arg

		 ti.d = pn + pv
	       unormuj(ti.d)
      jeżeli a2 < 0
		in = -n
	      a2 = iloczyn_skalarny(obsV, in)

	      arg = ((sinus(arcus_kosinus(a2)))*mD)/mU

	      jeżeli arg < 0.9999 to beta = arcus_sinus(arg)
      	wpp beta = 1.5705

	      arg = kosinus(beta)

		pn = n * arg
	      arg = mD/mU

		pv = (a2*in – obsV) * arg

  		ti.d = pn + pv
		unormuj(ti.d)
    ti.r = r.r + 1
		--mamy  ti kierunek załamania
   Jeżeli nie ma intersection(trójkąty, ti, nv, idx)
        Kolor_tła(ti, red, green, blue)
        jezeli r.c = RED to   fact2 = red
        (g,b,...)
	  rv = rv + nf * fact2
   wpp
	  InterpolujAlg2(trójkąty[idx], nv)
	  add = nf*rekurencyjny_kolor(ti, trojkąty, nv, nn, idx, nf)
	  rv = rv + add

	--kończenie obliczeń
 ret_val:
 Jeżeli rv > 1 to ZWRÓĆ 1
 ZWRÓĆ rv

Opisany powyżej algorytm jest uproszczony.


4. Opis struktur danych oraz ich zastosowania

4.1 Struktury danych.

Przedstawiam opis poszczególnych struktur danych.

Wierzchołek: LiczbaRzeczywista: x,y,z
	Używany jako punkt w przestrzeni 3D

Wektor: LiczbaRzeczywista: x,y,z
	Używany jako wektor w przestrzeni 3D, np. kierunek promienia, normalna itp.

Materiał: LiczbaRzeczywista: c,s,t
	Kolejne zmienne to: c – kolor, rozproszenie (color),
 s - odbijalność (specular),
 t – przezroczystość (transparency)
	Używany przy określaniu właściwości materiałów dla wierzchołków. Pomiędzy wierzchołkami interpolacja. Poza tym w wierzchołkach mogą też być koordynaty tekstur.

WspółczynnikTekstury: LiczbaRzeczywista: s,t
	Koordynaty tekstury: s,t e [0,1]

Powierzchnia: LiczbaRzeczywista: A,B,C,D
	Dokładny opis znajduje się w 2.4.4, przy opisie generowania. Używana przy obliczeniach przecięć, każdy trójkąt musi mieć obliczoną powierzchnię w której się zawiera.

Promień:
	Wierzchołek: P		- punkt od którego promień podąża
	Wektor: d			- kierunek w którym promień podąża
	LiczbaCałkowita: r,c 	- r to aktualny poziom rekursji, c to kolor (r lub g lub b)

	Promienie są wypuszczane w kierunku ekranu od obserwatora. Cały algorytm polega na dokładnym obliczeniu ich tras.


Trójkąt: Struktura została opisana dokładnie w 2.4.3, przy opisie trójkąta GPT.
	Całe sceny składają się z trójkątów, powierzchnie NURBS są triangulowane w preprocessingu i także ostatecznie są listą trójkątów.

Box: Axis Aligned Bounding Box (AABB), prostopadłościan wyrównany do osi układu, dokładny opis w 3.2. Drzewo boxów jest obliczane w preprocesingu lub wczytywane z pliku.
	Używany jako hierarchia otaczająca (BV – Bounding Volume)

Btree: Drzewo AABB, opis i sposoby generacji w 3.2
	Używane w celu przyspieszenia obliczeń raytacingu jako struktura lokalizująca kolizje.

Lista:
	WskaźnikNaListę: poprzedni, następny
	WskaźnikNaDane: dane

	Struktura używana do przechowywania list obiektów, które należy przekształcić, oraz bardzo często przechowująca dane tymczasowe różnych algorytmów.

Macierz: elementy[][]
	Przechowuje przekształcenia i transformacje świata, światła, obiektów, list obiektów itp. Alokowana dynamicznie, wielkość m i n należy określić alokując macierz.

ListaTransformacji:
	Macierz: m, mn
	LiczbaCałkowita index1, index2
	LiczbaCałkowita NURBS_index1, NURBS_index2
	Macierz m przechowuję transformację wierzchołków, a macierz mn transformację normalnych (pochodną). Indeksy 1 i 2 oznaczają od którego do którego trójkąta/NURBSa obowiązuje przekształcenie. Indeksy dla NURBS’a mogą być równe –1 gdy przekształcenie nie dotyczy NURBS’ów

Tekstura:
	LiczbaCałkowita długość, szerokość
Macierz: piksele,  każdemu pikselowi odpowiada kolor w formacie RBG, wielkość macierzy: długość x szerokość
Tekstury są używane do teksturowania obiektów na scenie, główny ekran także jest teksturą o wymiarach identycznych z wielkością okna.


NURBS: patrz też definicja NURBS w 2.1
	LiczbaCałkowita: idx – indeks(numer) powierzchni NURBS
	LiczbaCałkowita: ilośćTrójkątów – obliczane w trakcie triangulacji
	LiczbaCałkowita: tekstura (index tekstury na powierzchni NURBS)
	LiczbaCałkowita: n1, n2 (ilości punktów kontrolnych w kierunku u,v)
	LiczbaCałkowita: p1,p2 (wymiary w kierunku u,v)
	LiczbaCałkowita: m1,m2 (m = n + p + 1, wielkości tablic węzłów (knots) u i v)
	LiczbaCałkowita: d1,d2 (na ile trójkątów podzielić w kierunku u i v)
	TablicaLiczbRzeczywistych: t1,t2   węzły (nodes) t1[n1], t2[n2]
	TablicaLiczbRzeczywistych: knot1,knot2   węzły (knots) knot1[m1+1], knot2[m2+1]
	Macierz: w wagi poszczególnych punktów kontrolnych: w[n1][n2]
	Macierz: P Punkty kontrolne: P[n1][n2]

	Używana aby wczytać definicję powierzchni a następnie ją ztriangulować i umieścić na liście trójkątów

Woksel:
	Lista trójkąty	lista trójkątów w tym wokselu
	Box box		box otaczający woksel (dokładna wielkośc woksela)
	LiczbaCałkowita ilość_boxów ilo boxów ma woksel
	Lista woksele  woksele które są potomkami tego woksela

	4.2 Zastosowanie struktur danych w algorytmach

Światło jest przechowywane jako Wierzchołek (światło punktowe) lub jako Wektor: światło wektorowe.

Wektor/Wierzchołek light


Obraz generowany jest przechowywany jako tekstura o wielkości okna:

Tekstura obraz


Pozostałe tekstury są wczytywane z plików BMP/JPG i po wczytaniu wszystkich są przechowywane w tablicy tekstur:

Tekstura tekstury[]


Drzewo lokalizacyjne (AABB drzewo) jest trzymane w liście poszczególnych boxów. Wszystkie one tworzą poza tym strukturę drzewiastą. W liście mamy dostęp kolejno do wszystkich boxów, drzewo zapewnia szybką lokalizację. Te dwie struktury współdzielą dane.

ListaBoxów: boxy
AABBDrzewo: korzeń


Wszystkie trójkąty są przechowywane na liście, powierzchnie NURBS po ztriangulowaniu dodają „swoje” trójkąty do listy, dodają także indeks na nie wskazujący.

Lista trójkąty


Przekształcenia świata są trzymane w dwóch macierzach: jedna przechowuje globalne przekształcenie całego świata, druga globalne przekształcenie wszystkich normalnych.

Macierz M, MN


Przekształcenia list trójkątów są trzymane na liście. Jest to więc lista list trójkątów do przekształceń, szczegóły jak jest ona używana są w: 5.2.1.2

Lista listy_przekształceń_trójkątów

5. Zaimplementowane funkcjonalności, programy

	5.1 Odczyt i zapis danych

	Opiszę teraz obsługiwane przez mój algorytm formaty plików wejściowych i wyjściowych, a także które pliki do czego są używane.

5.1.1 Odczyt danych

	Pliki odczytywane przez mój algorytm można podzielić na 3 podstawowe typy:
	-pliki zawierające definicję sceny
	-pliki graficzne zawierające tekstury lub częściowo wygenerowane obrazy
	-pliki zawierające preprocesing sceny: zapisane AABB drzewo

	Opiszę poszczególne typy i sposoby ich wczytywania:

5.1.1.1 Wczytywanie definicji sceny

Program potrafi odczytywać pośrednio lub bezpośrednio sceny zapisane w następujących formatach:

DAT/FDAT: jest to standardowy tekstowy format pliku wejściowego używanego przez moją aplikację, może być wczytany bezpośrednio. Wszystkie inne formaty są konwertowane bezpośrednio do tego formatu. Format FDAT powstaje przez zapisanie bezpośrednio wszystkich trójkątów sceny po przekształceniach i triangulacji w oddzielnym pliku, używa się do tego opcji -B programu RAYS. Format pliku FDAT jest identyczny jak format pliku DAT, pliki te nie zawierają definicji powierzchni NURBS, tylko trójkąty składowe tych powierzchni po triangulacji. Oto przykładowy plik DAT zawierający wszystkie dostępne opcje wraz z komentarzami.

	Przykładowy plik z opisem wszelkich opcji jako komentarze znajduje się w pliku dat\options.DAT
	BIN: Format wewnętrzny mojego algorytmu w którym cała scena jest zapisana binarnie (bez przekształceń, po prostu każdy trójkąt po obliczeniu ostatecznej pozycji). Format ten może być wczytany bezpośrednio, pliki BIN można wygenerować za pomocą mojego algorytmu na podstawie plików DAT lub FDAT.

NURBS: pliki te zawierają definicję powierzchni NURBS, aktualna wersja RAYS
potrafi juz  odczytywać definicje NURBS w plikach DAT, więc pliki NURBS nie są już potrzebne - ich zawartość można bezpośrednio wstawiać do DAT. Aby odczytać plik NURBS należy go najpierw skonwertować na DAT za pomocą programu NURBS2DAT.

	IGES/IGS: Można odczytać definicje powierzchni NURBS z plików Iges, służy do tego program IGES2DAT konwertuje on plik IGES do pliku DAT, odnajduje w nim rekordy 126 i zamienia je na rekordy NURBS, można podać parametry dodatkowe powierzchni tj. gęstość triangulacji, kolor, skalowanie itp.

	TRI/ULI: pliki te przechowują trójkąty oraz normalne i koordynaty tekstur. Jest to format używany przez nas na VR i Grafice. Program ULI2DAT konwertuje pliki TRI/ULI na DAT, podobnie jak IGES2DAT można podać wiele opcji.

	MD2: format modeli w Quake2, brakuje w nim normalnych. Program MD22DAT konwertuje te pliki do DAT, także można podać wiele opcji.

	3DS: Pliki 3D Studio Max mogą być odczytane na dwa sposoby. Pierwszy program 3DS2TRI (starszy i gorszy) potrafi je skonwertować do formatu TRI/ULI, tracąc przy tym informacje o teksturach. Drugi (nowszy i lepszy) potrafi je skonwertować bezpośrednio do formatu DAT, zachowując informacje o teksturach. Można wybrać interpolację normalnych lub obliczenie ich zlecić mojemu algorytmowi opisanemu w 2.4.5. Program jest wysoce konfigurowalny, można zmienić wszystkie właściwości wierzchołka, zapisywać transformacje itp.

5.1.1.2 Wczytywanie tekstur i scen wygenerowanych częściowo

Program potrafi odczytywać tekstury w formacie BMP (Bitmapa 24bitowa). Jeżeli w czasie kompilacji ustawiono odpowiednią opcję to jest możliwość wczytywania formatu JPG/JPEG. Poniżej zakładam, że obsługa JPEG jest włączona.

Odczyt tekstury:

Tekstury są odczytywane według ich ID. Przeszukiwany jest katalog tekstur (który może być podany w definicji sceny lub z wiersza poleceń programu), najpierw następuje próba odczytu ID.bmp a potem (jeżeli odczyt BMP się nie powiódł) ID.jpg i ID.jpeg. Jeżeli ID jest równe 0 to algorytm zakłada, że obiekt nie ma tekstury. Odczyt plików JPEG odbywa się za pomocą napisanego przeze mnie algorytmu używającego biblioteki JPEGLIB, algorytm: load_jpeg_file w rayslib.c.

Odczyt sceny częściowo wygenerowanej:

Raytracing można przerwać w dowolnym momencie przez naciśnięcie CTRL+C, lub wysłanie sygnału z GUI, lub przez zdalne polecenie przerwania. Wysyłanie sygnałów i zdalne polecenia są opcjonalne. O tym czy są dostępne decydują flagi kompilacji. Program zapisuje wtedy to co obliczył do tej pory (z dokładnością do linii pionowej), oczywiście rozdzielczości obrazu wejściowego i obrazu generowanego muszą się zgadzać, ale jest to jedyna metoda sprawdzania czy kontynuujemy raytracing na tej samej scenie czy na innej. Jest więc możliwość stworzenia obrazu wynikowego składającego się z częściowych renderingów różnych scen.

Taki plik można odczytać i kontynuować obliczenia od miejsca gdzie je zakończono (z dokładnością do linii pionowej). Algorytm przechodzi przez kolejne pionowe linie obrazu aż znajdzie cała czarną linię – od tej linii zaczyna obliczenia. Ważne aby plik częściowo wygenerowany był bitmapą a nie JPEGiem. W pliku JPEG takie znaczniki jak cała linia czarna mogą zostać zniszczone przez kompresję/dekompresję stratną. W przypadku próby kontynuacji na podstawie pliku JPEG zostanie wyświetlone stosowne ostrzeżenie.

Algorytm wykrywania skąd należy kontynuować obliczenia:

ZnajdzOstatiąLinię(t)

Wejście: Tekstura tekstura
Wyjście: indeks linii od której zacząć obliczenia

dla i = 0 do eksturat.X wykonaj:
	ok = 1
	dla j = 0 do tekstura.Y wykonaj:
		jeżeli tekstura[i][j] <> czarny to ok. = 0
		ustaw_kolor(ekran, i, j, tekstura[i][j])
	jeżeli ok.
		ZWRÓC i

	5.1.1.3 Wczytywanie drzewa AABB

	Plik zawierający drzewo AABB ma rozszerzenie BTREE, plik ten może być binarny lub tekstowy (rozszerzenie się nie różni). Zalecany format to binarny ponieważ drzewa lokalizacji i tak nie edytuje  się „ręcznie”, plik binarny jest szybciej wczytywany. Procedura zapisu i odczytu takiego pliku została opisana w 3.2.2.1 (zapis) i 3.2.2.2 (odczyt)

	5.1.2. Zapis danych

Ogólnie zapis danych możemy podzielić na 5 grup:

	-zapis obrazów wynikowych/częściowych
	-zapis binarnych scen
	-zapis preprocesingu
	-zapis przekształcenia świata/światła w trybie podglądu
	-konwersja jednego formatu wejściowego na drugi

5.1.2.1 Zapis obrazów wynikowych/częściowych

	Obrazy wynikowe są zapisywane w formacie BMP a także jeżeli odpowiednio skonfigurowano program to w formacie JPEG. Istnieje także możliwość zapisania JPEG’ów w odcieniach szarości a także podwojonych rozmiarów bitmap w przypadku użycia antyaliasingu (potrzebne by potem wznowić obliczenia). Obrazy są zapisywane w różnych sytuacjach, takich jak:

	-zakończenie obliczeń
	-błąd krytyczny
	-sygnał przerwania (jeżeli wkompilowano obsługę sygnałow)
	-żądanie użytkownika
	-co pewną (pokreśloną przez użytkownika) ilość linii

	Dla każdej z tych sytuacji tworzone są inne pliki BMP/JPEG. Jest możliwość konfigurowania nazw tych plików a także tego kiedy będą generowane (np. co ile linii – autobackup). Opcję autobackup można skonfigurować w pliku ze sceną lub z wiersza poleceń programu, a także przez komunikację z serwerem przez internet lub przez odpowiedni klawisz na oknie GUI.. Obraz częściowy jest zapisywany do poprzedniej linii w stosunku do aktualnie obliczanej, aktualna jest zapisywana w całości na czarno – jako znacznik końca danych obliczonych, patrz odczyt częściowego obrazu w 5.1.1.2. Algorytm zapisu wygląda następująco:

WyczyśćLinie(tekstura,idx)

Wejście: tekstura, idx: indeks aktualnej linii pionowej
Wyjście: tekstura przygotowana do zapisu częściowego

Od i = 0 do tekstura.Y wykonaj:
	Tekstura[idx][i] = CZARNY


5.1.2.2 Zapis binarnych scen

	Sceny binarne BIN, mogą być wygenerowane z plików wejściowych DAT, jeżeli użytkownik poda odpowiednią opcję z wiersza poleceń programu. Pliki BIN mogą być bezpośrednio wczytywane przez algorytm, patrz: 5.1.1.1. Nie jest to jednak format zalecany, gdyż nie istnieje praktycznie żadna możliwość ingerencji w scenę, format binarny jest zalecany dla drzew AABB, których ręczna edycja na ogół nie ma sensu.


5.1.2.3 Zapis obliczonego drzewa AABB

	Wygenerowane drzewo AABB może być zapisane w pliku BTREE, jeżeli zostanie podana jedna z dwóch opcji zapisywania drzewa z wiersza poleceń aplikacji. Jedna z opcji nakazuje zapis binarny drzewa a druga zapis tekstowy. Tylko jedna opcja na raz jest dozwolona, w przypadku użycia obu, obowiązującą będzie ta podana na końcu. Zapis drzewa BTREE jest opisany w 3.2.2.1 a odczyt w 3.2.2.2 i 5.1.1.3

5.1.2.4 Zapis przekształceń w trybie podglądu

	W trybie podglądu (renderingu normalnego a nie raytracingu), można dowolnie przekształcać scenę i światło. Po naciśnięciu odpowiedniego przycisku aktualne przekształcenie sceny i światła zostaje zapisane w pliku world_trans.DAT. Nie jest to poprawny plik DAT (nie zawiera nagłówka i danych), ale przekształcenie w nim zapisane można zapisać w pliku sceny na którym działamy za pomocą zwykłego edytora tekstu np. VI.. Przekształcenie jest zapisywane jako ListTransform/WorldTransform i LightTransform.
	5.1.2.5 Konwersje między formatowe

	Dodatkowe programy narzędziowe (poza RAYSLib) są w stanie konwertować pewne formaty danych do innych wczytywalnych przez RAYSLib: DAT, FDAT, BIN, BTREE. Krótkie omówienie konwersji danych jest w 5.1.1.1. Dokładniejsze informacje na ten temat można znaleźć w rozdziale 6, przy opisie poszczególnych programów.

5.2 Przekształcenia sceny i kopiowanie obiektów

	W programie istnieje możliwość przekształcania świata, światła, trójkątów, list trójkątów, NURBS’ów. Istnieje także możliwość kopiowania trójkątów, list trójkątów, NURBS’ów, list NURBS’ów. Przy kopiowaniu obiektów można cofnąć zadane dla nich transformacje, a także grupowo zadać pewne nowe właściwości: np. nową teksturę.

5.2.1 Przekształcenia obiektów

	Przekształcenia obiektów możemy podzielić na:
		-przekształcenia świata
		-przekształcenia światła
		-przekształcenia obserwatora/kamery
		-przekształcenia pojedynczych trójkątów
		-przekształcenia list trójkątów
		-przekształcenia powierzchni NURBS
		-przekształcenia listy NURBS’ów

	Kolejność stosowania przekształceń. Najpierw konkretne dla danego trójkąta (jeżeli są), potem przekształcenia listy a na końcu przekształcenia świata. Przekształcenia NURBS’ów są niejawnie konwertowane do przekształceń odpowiednich list trójkątów, tj tych trójkątów które do danej powierzchni (lub lity powierzchni) należą. Przekształcenia kamery i światła są niezależne od opisanych powyżej przekształceń. Jeżeli jest kilka przekształceń list trójkątów to: przekształcenia list trójkątów są stosowane w kolejności od ostatniego zadanego (najniżej w pliku) do pierwszego (najwyżej w pliku)

5.2.1.1 Przekształcenia pojedynczego trójkąta

W celu przekształcenia trójkąta musimy dysponować macierzą przekształcenia wierzchołków i normalnych (pochodna macierzy przekształcenia wierzchołków). Algorytm przekształcający pojedynczy trójkąt jest następujący

PrzekształćTrójkąt(t, m, mn)

Wejście: trójkąt, macierze przekształceń wierzchołków i normalnych
Wyjście: przekształcony trójkąt

T.a = T.a * m
T.b = T.b * m
T.c = T.c * m

T.na = T.na * mn
T.nb = T.nb * mn
T.nc = T.nc * mn

	Przedstawiam jeszcze fragment algorytmu obliczającego macierze przekształceń dla translacji, rotacji, skalowania oraz podawanie własnej macierzy. Algorytm ten analizuje przekształcenie i odpowiednio oblicza macierz, w rzeczywistości obsługiwane jest dużo więcej przekształceń, patrz read_transformation_long w rayslib.c


OdczytajTransformacje(plik, m, mn)

Wejście: Plik z którego czytamy, macierze przekształceń wierzchołków i normalnych
Wyjście: Obliczone macierze przekształceń wierzchołków i normalnych

odczytaj transformacje z pliku
Jak odczytano „Translate”
	odczytaj x,y,z
	Translacja(m, x,y,z)
Jak odczytano „Scale”
	odczytaj x,y,z
	Skalowanie(m, x,y,z)
	Skalowanie(mn, x,y,z)
Jak odczytano „Rotate”
	odczytaj x,y,z
	RotacjaX(m,  x)
	RotacjaX(mn, x)
	RotacjaY(m,  y)
	RotacjaY(mn, y)
	RotacjaZ(m,  z)
	RotacjaZ(mn, z)
Jak odczytano „Macierz”
Odczytaj m, mn (tutaj można podać zupełnie dowolne przekształcenie wierzchołków i normalnych)
KONIEC

Translacja(m, x,y,z)

Wejście: Macierz: m, translacja: x,y,z
Wyjście: obliczona macierz m

M = I
M[0][3] = x
M[1][3] = y
M[2][3] = z

Skalowanie(m, x,y,z)

Wejście: Macierz: m, skalowanie: x,y,z
Wyjście: obliczona macierz m

M = I
M[0][0] = x
M[1][1] = y
M[2][2] = z

RotacjaX(m, kąt)

Wejście: Macierz: m, kąt obrotu
Wyjście: obliczona macierz m

 m[1][1] = cos(kąt)
 m[2][1] = sin(kąt)
 m[1][2] = -sin(kąt)
 m[2][2] = cos(kąt)

RotacjaY(m, kąt)

Wejście: Macierz: m, kąt obrotu
Wyjście: obliczona macierz m

 m[0][0] = cos(kąt)
 m[2][0] =- sin(kąt)
 m[0][2] = sin(kąt)
 m[2][2] = cos(kąt)

RotacjaZ(m, kąt)

Wejście: Macierz: m, kąt obrotu
Wyjście: obliczona macierz m

 m[1][1] = cos(kąt)
 m[0][1] = sin(kąt)
 m[1][0] = -sin(kąt)
 m[0][0] = cos(kąt)

Przekształcenia kamery, światła i świata są analogiczne, z tym że kamera nie ma normalnej a światło jeżeli jest punktowe to jest przekształcane przez M a jeżeli wektorowe to przez MN. Przekształcenia świata zmieniają globalną macierz świata, która jest stosowana w ostatniej kolejności dla wszystkich obiektów na scenie.

5.2.1.2 Przekształcenia list trójkątów.

	Na liście do przekształceń znajdują się macierze przekształceń oraz indeksy: dolny i górny trójkątów, które będą przekształcane. Następujący algorytm przegląda wszystkie listy przekształceń i sprawdza dla każdej listy czy trójkąt do niej należy. Jeżeli należy, to stosowane jest dla niego przekształcenie tej listy. Wynika z tego, że przekształcenia list trójkątów są stosowane w kolejności od ostatniego zadanego (najniżej w pliku) do pierwszego (najwyżej w pliku). Aby przekształcić wszystkie trójkąty należy dla każdego z nich wywołać następujący algorytm. Ostateczne przekształcenie każdego trójkąta, który musi być przekształcony realizuje funkcja opisana w 5.2.1.1

PrzekształćTrójkątListamiPrzekształceń(t)

Wejście: trójkąt przekształcany
Wyjście: przekształcony trójkąt

Tmp = głowa_listy_list_przekształcających

Dopóki są listy przekształceń w „tmp”

	Lt  = pobierz_listę_przekształceń(tmp)
	Jeżeli t.idx >- lt.indeks_dolny i t.idx <= lt.indeks_górny
		PrzekształćTrójkąt(t, lt.M, lt.MN)
	Tmp.następny


Dokładana lista możliwych przekształceń znajduje się w pliku z przykładami użycia wszystkich opcji: dat/options.DAT. Przekształcenia możemy podzielić na:


	-afiniczne: 		opisane macierzami wierzchołków i normalnych

-materiałowe:	przekształcenia zmieniające właściwości materiału, takie jak: tekstura, kolor, odbijalność, przepuszczalność, rozbłysk itp. Dla tych przekształceń można zadać warunek, że tylko trójkąty z zadaną tekstura podlegają przekształceniu

5.2.1.3 Przekształcenia list powierzchni NURBS

	Najpierw dla zadanych indeksów NURBS’ów obliczane są trójkąty w których są zawarte te powierzchni. Dla tych trójkątów tworzona jest lista przekształcenia jak opisałem w 5.2.1.2. Oto algorytm obliczający indeksy początkowego i końcowego trójkąta:


ObliczIndeksyTrójkątów(i, j)

Wejście: indeksy powierzchni NURBS dolny i górny w i, j
Wyjście: indeksy trójkątów dolny i górny w i,j

 i1 = i
 i2 = j

 Jeżeli i1 lub i2 mniej niż zero to BŁĄD
 Jeżeli i1 lub i2 więcej niż ilość NURBS’ów to BŁĄD
 Jeżeli j < i to BŁAD

 i = nurbs[i1].idx
 j = nurbs[i2].idx + nurbs[i2].ntri –1



	Obliczanie indeksów trójkątów powierzchni NURBS w preprocesingu

Indeksy trójkątów w powierzchni NURBS są przy ich triangulacji, trójkąty powierzchni NURBS są dodawane za trójkątami sceny, oto algorytm obliczający indeksy trójkątów powierzchni NURBS, każda powierzchnia NURBS ma już obliczone ile ma trójkątów (policzono to w trakcie triangulacji: d1 x d2 x 2)

ObliczanieIndeksówNURBSów(nNURBS, nTris)

Wejście: ilość NURBSów i trójkątów
Wyjście: wszystkie powierzchnie NURBS z obliczonymi indeksami trójkątów
	    Całkowita ilość trójkątów w powierzchniach NURBS

Nt = 0

Dla i = 0 do nNURBS wykonaj:
	Nurbs[i].idx = nT + nTris
	Nt = nt + nurbs[I].nTris

ZWRÓC nt

Przy tak obliczonych indeksach możemy przekształcać powierzchnie (lub listy powierzchni) NURBS jak zwykłe listy trójkątów co zostało opisane w 5.2.1.2

5.2.2 Kopiowanie obiektów

Mój algorytm zakłada także, że możliwe jest kopiowanie obiektów. Ogólnie kopiowanie obiektów możemy podzielić ze względu na typy i ilość obiektów kopiowanych na:

	-kopiowanie pojedynczego trójkąta
	-kopiowanie listy trójkątów
		-z wycofywaniem różnych przekształceń
		-bez wycofywania przekształceń
	-kopiowanie pojedynczej powierzchni NURBS
	-kopiowanie wielu powierzchni NURBS

	Przykłady różnych kopiowań obiektów znajdują się w pliku dat/options.DAT. Omówię teraz poszczególne rodzaje kopiowania i algorytmy użyte aby je zaimplementować..


	5.2.2.1 Kopiowanie pojedynczego trójkąta

	Najpierw sprawdzane są indeksy źródła i docelowy aby zweryfikować czy są poprawne, w zależności od tego czy wyspecyfikowano cofnięcie przekształceń czy nie jest ono wykonywane.

	Najczęściej kopiujemy grupę trójkątów tworzącą np. jedna ścianę jakiejś bryły. Powielamy ją tyle razy ile jest ścian. Aby efekt był dobry należy wycofać np. przekształcenia świata aby łatwiej przenieść tą ścianę w odpowiednie miejsce. Często też kopiujemy już przekształcone listy trójkątów aby dokonać na nich innego przekształcenia, wtedy należy obowiązujące przekształcenie wycofać. Dlatego dodałem możliwość wybrania tego czy cofamy przekształcenia świata i listy przekształceń czy nie. Po cofnięciu listy przekształceń można zastosować inną na nowo powstałej grupie skopiowanych trójkątów. Przekształcenia własne danego trójkąta nie są cofane, gdyż należą one tylko do niego i są dla niego typowe (nawet nie są nigdzie zapamiętywane, tylko obliczane raz na stałe przy tworzeniu trójkąta – są zapisane na stałe w jego współczynnikach). Przedstawię teraz algorytm kopiujący pojedynczy trójkąt:

KopiujTrójkąt(z, do)

Wejście: indeksy trójkątów źródła i docelowego
Wyjście: skopiowany trójkąt

Jeżeli do = z to BŁĄD
Jeżeli z > do to BŁĄD (można kopiować tylko z już odczytanych)
Jeżeli z <0 lub z > ilość trójkątów to BŁĄD
Trójkąt[do] = Trójkąt[z]
Trójkąt[do].idx = z
Jeżeli były przekształcenia i użytkownik chce to wycofaj_przekształcenia_świata(Trójkąt[do])
Jeżeli były przekształcenia i użytkownik chce to wycofaj_listę_przekształceń(Trójkąt[do])
PrzekształćTrójkątListamiPrzekształceń(Trójkąt[do])
PrzekształćTrójkąt(Trójkąt[do], macierzeŚwiata M i N)

5.2.2.2 Kopiowanie listy trójkątów

	Przekształcenie listy trójkątów jest podobne do przekształcenia pojedynczego trójkąta: na początku sprawdzane są indeksy kopiujące, a potem jest wywoływana metoda dla pojedynczego trójkąta. W metodzie kopiujące wiele trójkątów na raz można określić czy wycofywać przekształcenia świata i list czy nie. Oto algorytm:

KopiujTrójkąty(do, z, ile, wycofajListe, wycofajŚwiat)

Wejście: indeksy do, z, ilość, oraz czy wycofywać poszcz. operacje
Wyjście: skopiowane „ile” trójkątów

Jeżeli ile < 0 to BŁĄD
Jeżeli src + (ile-1) >= dst to BŁAD (Wejdziemy na aktualny lub jeszcze nie wczytany)
Jeżeli src < 0 lub src > ilośćTrójkątów-(ile-1) to BŁĄD
Dla i = do dopóki  i < do+ile wykonaj:
	KopiujTrójkąt(src+I-do, I)

	5.2.2.3 Kopiowanie powierzchni NURBS

Metoda kopiująca powierzchnię NURBS jest wywoływana w preprocessingu. Kopiuje ona cała strukturę NURBS do innej docelowej (sprawdziwszy najpierw indeksy tak jak w 5.2.2.1). procedura jest wywoływana przed triangulacją, więc oddzielne trójkąty powstaną dla skopiowanej i oryginalnej powierzchni NURBS. Potem można obie powierzchnie przekształcać zarówno ListTransform (czyli podając ich trójkąty, patrz: 5.2.1.2) jak i NURBSTransform (czyli całe powierzchnie NURBS, patrz 5.2.1.3). Metoda kopiująca nie wycofuje żadnych przekształceń zadanych na powierzchni NURBS, ale można swobodnie stosować przekształcenia zarówno dla powierzchni źródłowej jak i docelowej, ponieważ kopiowanie odbywa się przed jakimkolwiek użyciem transformacji, nawet przed triangulacją!.
Nie istnieje metoda kopiująca NURBS’y grupowo, pewnym zastępstwem tej metody może być kopiowanie listy trójkątów zawierającej zadane powierzchnie NURBS metodą CopyTriangles lub CopyTrianglesAdv, opisanej w 5.2.2.2. Trzeba jednak znać indeksy początkowe i końcowe tej listy – w czasie triangulacji wyświetlane są indeksy i ilości trójkątów należących do poszczególnych powierzchni NURBS.


5.3 Inne algorytmy

5.3.1 Obsługa komentarzy w pliku wejściowym

	W pliku wejściowym DAT mogą się znajdować komentarze w stylu C tj:
/* komentarz */ i // komentarz. Pierwszy zaczyna się od /* a kończy na */, a drugi zaczyna się od // i kończy wraz z końcem linii.

5.3.2 Obsługa sterowania przez Internet

	Zaimplementowany został prosty serwer nasłuchujący polecenia na zadanym porcie. Za pomocą dowolnego klienta (np. telnet), można połączyć się z serwerem i wydawać polecenia. Dokładny opis działania serwera i poleceń znajduje się w rozdziale 6. Między innymi można zmieniać głębokość rekursji, generować tymczasowe obrazy, zmieniać chropowatości powierzchni, wyświetlać statystyki raytracingu i wiele innych.. Serwer działa na oddzielnym wątku, którego jedynym zadaniem jest nasłuch połączeń. Serwer jest prosty i nie obsługuje wielu klientów na raz, jest to serwer iteracyjny. To czy serwer jest wkompilowany czy nie zależy od tego jakie flagi kompilacji zostały zastosowane.


5.3.3 OpenGL GUI

	Jest możliwość wkompilowania w algorytm interfejsu graficznego w OpenGL. Interfejs graficzny działa w oddzielnym wątku i służy tylko do wyświetlania efektu częściowego oraz przyjmowania poleceń. Jeżeli interfejs ten został wkompilowany to możemy uruchomić program w dwóch trybach:

	-tryb raytracingu z wyświetlaniem w OpenGL
	-tryb podglądu sceny w OpenGL

	5.3.3.1 Wyświetlanie Raytracingu w OpenGL

	Przeprowadzany jest normalny raytracing, wynik można oglądać na bieżąco w oknie graficznym OpenGL. Częstotliwość odświeżania okna może być konfigurowana przez użytkownika, w oknie można też naciskać różne przyciski wysyłając w ten sposób polecenia do algorytmu raytracingu. Wykaz poleceń i ich dokładne działanie jest opisany w rozdziale 6. Między innymi można zmieniać głębokość rekursji, generować tymczasowe obrazy, zmieniać chropowatości powierzchni, wyświetlać statystyki raytracingu i wiele innych

	5.3.3.2 Podgląd sceny w OpenGL

	W tym trybie scena nie jest obliczana przez raytracing, ale jej uproszczenie jest renderowane przez OpenGL, zysk jest taki, że możemy w czasie rzeczywistym przekształcać scenę: translacje, rotacje, skalowania. Możemy także włączać/wyłączać: światło, tekstury, przezroczystość itp. Możemy oglądać punkty kontrolne NURBS’ów, zmieniać tryb renderowania z wireframe na solid i odwrotnie. Opis wszystkich opcji podglądu znajduje się w rozdziale 6. Najważniejszą opcją jest możliwość zapisu przekształceń w pliku world_trans.DAT i potem użycie go przy raytracingu. Patrz opis w 5.1.2.4

5.3.4 Chropowatość powierzchni

	Istnieje możliwość zadania chropowatości wybranym (lub wszystkim) trójkątom. Dodano w tym celu właściwość NormalDistorber (RozproszenieNormalnej). Jest to właściwość każdego trójkąta. Aby rozpraszanie normalnych było brane pod uwagę należy uruchomić program z odpowiednią opcją, opis znajduje się w rozdziale 6. Kolejność obliczania współczynnika rozproszenia dla trójkąta. Najpierw jest brane rozproszenie podane w WorldTransform (jeżeli jest). Jeżeli podano globalne  rozproszenie: opcja z wiersza poleceń użytkownika (dokładny opis w rozdziale 6), lub globalne w pliku to zastępuje ono to z WorldTransform.. Jeżeli trójkąt ma własne to go używa, a jeżeli jest przekształcenie listy, to zastępuje ono to zdefiniowane dla trójkąta.

Algorytm obliczający rozproszenie jest prosty. Po obliczeniu normalnej zostanie ona zaburzona o losową wartość z przedziału [0,współczynnik_zaburzenia]. Zaburzane są także koordynaty tekstury (materiał jest chropowaty).Oto algorytm:

RozproszNormalną(wsp, n, tc)

Wejście: wsp: wspólczynnik rozproszenia, n – normalna, tc – koordynaty tekstury
Wyjście: zburzona normalna i koordynaty tekstury

Dist = rand(wsp)
Dist = 2*dist – 1

n = n + dist

tc = tc + f(dist)	(jakaś funkcja zaburzenia, np f = 0.2)

Jeżeli tc < 0 to tc = 0
Jeżeli tc > 1 to tc = 1

Unormuj(n)

	5.3.5 Obsługa sygnału błędu segmentacji.

	Sygnał błędu segmentacji (SIGSEGV), jest przechwytywany i wersja debug programu wchodzi w nieskończoną pętle, oczekując na przyłączenie debuggera. Można wtedy sprawdzić CallStack i łatwiej znaleźć przyczynę wystąpienia błędu segmentacji (naruszenia ochrony pamięci).

	5.3.6. Generowanie sceny w odcieniach szarości.

	Można ustawić algorytm aby obliczał kolory tylko dla jednego promienia (zielonego) i wszystkim pozostałym przypisywał to samo natężenie R=G=B. Wtedy otrzymamy scenę w odcieniach szarości, ale 3x szybciej. Wszelkie efekty rozszczepienia światła (pryzmaty itp.) będą niewidoczne.
	5.3.7 Łączenie i transformacje AABB-drzew

Program BTREECONV obsługuje transformowanie i łączenie AABB drzew. Przedstawię używane algorytmy:

5.3.7.1 Algorytm łączenia AABB-drzew:

PołączAABBdrzewa(t1, t2)
Wejście: AABBDrzewo t1,t2 (najlepiej rozłączne)
Wyjście: AABBDrzewo t połączone drzewa t1 i t2

WcztajDrzewo(t1)
WczytajDrzewo(t2)
lista_trójkątów = połącz_listy(lista1, lista2)
dla i = 0 do ilość_trójkątów2 wykonaj:
	trójkąt[i].indeks = i + ilość_trójkątów1
t.korzeń.lewa_gałąź    = t1
t.korzeń.prawa_gałąź = t2
b = t.korzeń.box
b.minX = min(b1.b.minX, b2.b.minX)
b.minY = min(b1.b.minY, b2.b.minY)
b.minZ = min(b1.b.minZ, b2.b.minZ)
b.maxX = max(b1.b.maxX, b2.b.maxX)
b.maxY = max(b1.b.maxY, b2.b.maxY)
b.maxZ = max(b1.b.maxZ, b2.b.maxZ)
ZapiszDrzewo(t)

5.3.7.2 Algorytm przekształcania AABB-drzew:

TransformujDrzewo(t)
Wejście: t-AABB drzewo
Wyjście: przekształcone AABB drzewo

jeżeli ma lewą gałąź to TransformujDrzewo(t.lewa_gałąź)
jeżeli ma prawą gałąź to TransformujDrzewo(t.prawa_gałąź)
Box b = t.box
b.minX = b.minX * skalaX
(... pozostałe elementy boxu)
RotacjaBoxu(b)
b.minX = b.minX + translacjaX
(... pozostałe elementy boxu)


RotacjaBoxu(b)
Wejście: Box b
Wyjście: obrócony box b
(niech y oznacza miny i maxy, analogicznie x,z)
(dla klarowności pomijam fakt, że do zamiany powinna być użyta zmienna pomocnicza)


Jeżeli obrót wokół OX o 90   stopni to: y = z	z = -y
Jeżeli obrót wokół OX o 180 stopni to: y = -y	z = -z
Jeżeli obrót wokół OX o 270 stopni to: z = y	y = -z
Jeżeli obrót wokół OYo 90    stopni to: z = x	x = -z
Jeżeli obrót wokół OY o 180 stopni to: z = -z	x = -x
Jeżeli obrót wokół OY o 270 stopni to: x = z	z = -x
Jeżeli obrót wokół OZ o 90   stopni to: x = y	y = -x
Jeżeli obrót wokół OZ o 180 stopni to: x = -x	y = -y
Jeżeli obrót wokół OZ o 270 stopni to: y = x	yx= -y

jeżeli minx > maxx to minx <-> maxx
jeżeli miny > maxy to miny <-> maxy
jeżeli minz > maxz to minz <-> maxz


	1.	 Opis użycia wszystkich programów, dane techniczne

6.1 Dane techniczne

	Program jest napisany w języku C. Nic poza standardową biblioteką C nie jest wymagane, więc program można skompilować na praktycznie dowolnym systemie operacyjnym (najuboższą wersję). W zależności od dostępnych bibliotek można dołączyć biblioteki:

-LibJPEG: 		odczyt/zapis w formacie JPEG (obrazy i tekstury)
-POSIX Signals 	obsługa przerywania/wznawiania RT, zapisu na żądanie itd.
-POSIX Inet 	obsługa programu przez sieć
-OpenGL 		wyświetlanie w trakcie raytracingu

	Platforma macierzysta wszystkich aplikacji to FreeBSD UNIX 5.2, 5.3, 5.4, 6.0, 6.1
Procesor i386 i amd64.

	6.2 Opisy programów

Moja aplikacja raytracingu składa się z głównego programu realizującego algorytm raytracingu oraz szeregu programów narzędziowych służących do generacji sceny, konwersji między formatami oraz służących do innych celów które zostaną opisane poniżej. Przedstawię teraz opis głównej aplikacji oraz programów narzędziowych:

Programy składowe:

RAYS 		główny program do RayTracing'u

NURBS	program generujący powierzchnie NURBS za pomocą interpolacji (generuje powierzchnię zadanych stopni u i v przechodzącą przez
zadane punkty – metody numeryczne; interpolacja globalna)

NURBS2DAT	Zamienia pliku .NURBS na .DAT (aktualna wersja RAYS potrafi juz 				to robić automatycznie, więc pliki NURBS nie są już potrzebne - 					ich zawartość można bezpośrednio wstawiać do DAT)

IGES2DAT		konwertuje plik IGES do pliku DAT, odnajduje w nim rekordy 126 i 				zamienia je na rekordy NURBS, można podać parametry dodatkowe 				powierzchni tj. gęstość triangulacji, kolor, skalowanie itp.

BTREECONV	program służący do manipulacjami drzew AABB. Potrafi konwertować 				zapisane drzewa AABB z formatu binarnego do tekstowego, 					dokonywać przekształceń na drzewach (skalowania, translacje i 					ograniczone rotacje) a także potrafi łączyć drzewa.

ULI2DAT		konwertuje pliki TRI/ULI (używane przez nas na VR) na DAT, 					podobnie jak IGES2DAT można podać wiele opcji

MD22DAT		konwertuje pliki danych Quake'a (MD2) do DAT, można podać wiele 				opcji

3DS2TRI		konwertuje pliki 3DS do plików TRI (potem można użyć ULI2DAT by 				uzyskać plik DAT), normalne są obliczane poprzez interpolacje w 				wierzchołkach, ponieważ w plikach 3DS brak normalnych

3DS2DAT	konwertuje pliki 3DS do formatu DAT, możliwość interpolowania normalnych lub pozostawienia ich obliczenie dla RAYSLIB (powstaną płaskie trójkąty), zapisuje także przekształcenia i tekstury.

60FACES		generuje 60 ścian, podanie odpowiednich opcji generuje 						ztriangulowany 12-ścian

TORRUSGEN, CONE, CUBE, CYLINDER, BALL
generują odpowiednio: torus, stożek, sześcian, walec, kulę, można określić wiele parametrów brył

RANDNURB,RTRIANGLE
generują losowe powierzchnie NURBS i zbiory trójkątów

TABLE		generuje stół z teksturami, przezroczystością

TERMINAL, GETBMP:
służą do komunikacji z serwerem RAYS pierwszy umożliwia 		wysyłanie zapytań przez sieć do serwera, drugi pobiera obliczany 		obraz (jako BMP, o co dotychczas wygenerowano)

WRAPPER		interaktywny program pytający o wszystkie opcje i ostatecznie 					wywołujący odpowiednio RAYS

TEX			generuje przykładową teksturę w formacie BMP






6.2.1 Program RAYS

	Rays jest główną aplikacją raytracingu, w zależności od systemu operacyjnego na którym działamy nazwa programu to: cyg_rays.EXE (Windows) lub rays (UNIX/Linux). Przedstawiam poniżej opis działania i dostępne opcje programu rays.

6.2.1.1 Opis opcji wiersza poleceń

-i „nazwa pliku”	podanie nazwy pliku sceny, np. –i dat/options.dat, domyślnie „scene.dat”

-o „nazwa pliku”	podanie nazwy pliku wynikowego, np. –o wynik.bmp
			jeżeli użyjemy jeszcze opcji –J, -K, -g, to powstaną pliki:
			wynik.jpeg, wynik_gs.jpeg, wynik_aa.bmp, domyślnie „screen.bmp”

-U „nazwa pliku”	podanie nazwy pliku w którym zostaną zapisane informacje
			debugowe, np. –U debug.txt, domyślnie „debug.out”

-D „nazwa pliku”	nazwa pliku gdzie będą zapisywane wyniki pośrednie na życzenie
			użytkownika (reakcja na SIGUSR1), np. –D wynik.bmp
			jeżeli użyjemy jeszcze opcji –J, -K, -g, to powstaną pliki:
wynik.jpeg, wynik_gs.jpeg, wynik_aa.bmp, domyślnie „ondemand.bmp”

-S „nazwa pliku”	nazwa pliku gdzie będzie zapisany wynik pośredni w reakcji na przerwanie (CTRL+C lub polecenie użytkownika), reakcja
	na sygnał SIGINT, np. –S przerwane.bmp, domyślnie „signalled.bmp”

-P „nazwa pliku”	nazwa pliku gdzie zostanie zapisany dotychczasowy wynik gdy wystąpi błąd krytyczny (tzw. Panic) –P error.bmp, domyślnie “panic.bmp”

-p “nazwa pliku”	nazwa pliku gdzie zapisywane będą wyniki pośrednie co pewien ustalony przez użytkownika czas (dokładniej co ustaloną ilość linii), np. –p tymczasowa.bmp, domyślnie „partial.bmp”

-T „ścieżka do katalogu”
ustaw katalog w którym system będzie poszukiwał tekstur, nie jest to docelowy katalog z teksturami, a tylko miejsce z którego system będzie szukał „TextureDirectory”, domyślnie „.”, podanie np. /root/ będzie oznaczało, że w katalogu /root/ będzie szukany zadany w pliku sceny katalog np. texture, więc ostatecznie teksturty powinny być w /root/texture. Katalog powinien być zakończony „/”, np. –T /home/bla

-r liczba		maksymalny poziom rekursji, np. –r 5, domyślnie: 6

-R „nazwa pliku”	nazwa pliku z którego system ma spróbować odczytać już obliczone dane i kontynuować raytracing. Powinna to być bitmapa, np. –R tyczasowa.bmp, domyślnie ten argument nie jest używany.

-b liczba	określa co ile linii zapisywać plik tymczasowy (patrz –p), np. –p 256, domyślnie 64

-x liczba	ustawia rozdzielczość poziomą obrazu wynikowego, np. –x 1024, domyślnie używana jest opcja odczytana z pliku ze sceną

-y liczba	ustawia rozdzielczość pionową obrazu wynikowego, np. –y 768, domyślnie używana jest opcja odczytana z pliku ze sceną

-s liczba	ustaw „seed” randomu, domyślnie używany jest aktualny czas, np. –s 1742389779

-n procent	ustawia globalne zaburzenie normalnych (chropowatość materiału) na zadany procent, aby zaburzenia normalnych działały musi być użyta także opcja –N, np. –n 1.2, domyślnie 0

-t milisekund	ustawia czas po jakim następuje odświeżenie okna OpenGL, np. –t 1000 (co sekundę), domyślnie 500. Podanie zbyt małej wartości spowoduje duże obciążenie procesowa obliczaniem grafiki, a maksymalne FPS, może spowodować i tak większe opóźnienia.

-m liczba 	minimalny cień rzucany przez obiekt, liczba z przedziału [0,1]. ustawienie na –m 0.1 spowoduje, że cień 0.1 będzie rzucany nawet przez 100% przezroczyste obiekty, domyślnie 0.05

-M liczba	maksymalny cień rzucany przez obiekt, liczba z przedziału [0,1]
	ustawienie np. –M 0.75 spowoduje, że nawet idealnie nieprzezroczysty obiekt będzie rzucał cień 0.75, domyślnie 0.5

-a liczba	światło tła, liczba z przedziału [0,1], minimalne oświetlenie obiektu, światło rozproszone dochodzące ze wszystkich stron jednakowo, np. –a 0.35, domyślnie 0.3

-k liczba całkowita	decyduje o sposobie liczenia efektu Fresnela. Liczby 0,1,2 są specjalne: 0 – oznacza nie licz w ogóle, 1 – oznacza efekt liniowy (dobre przybliżenie), 2 – oznacza efekt kwadratowy (bardzo intensywny efekt). Pozostałe liczby są dzielone przez 1000 i używane jako wykładnik proporcjonalności, np. –k 1500 da proporcjonalność z potęgą 1.5, domyślnie 775.

-q procent	określa jakość JPEG’a przy zapisie za pomocą LibJPEG, np. –q 60, domyślnie 90.

-Q kolor RGB	kolor tła, format jest następujący RRGGBB, gdzie RR,GG,BB to 2-cyfrowe liczby hexadecymalne od 00 do FF określające nasycenie składowych czerwonej, zielonej i niebieskiej, np. –Q FFFF00 da nam żółty kolor, domyślnie: 7F7F7F: szary.

-W kolor RGB	kolor światła, format jest następujący RRGGBB, gdzie RR,GG,BB to 2-cyfrowe liczby hexadecymalne od 00 do FF określające nasycenie składowych czerwonej, zielonej i niebieskiej, np. –Q 0000FF da nam niebieski kolor, domyślnie: FFFFFF: biały.

-w	ustawia tryb tła (gdy nie jest używana tekstura tła). Ustawienie na 1 powoduje używanie koloru zadany przez –Q lub domyślnego. Dla –w 0 kolor domyślny lub zadany przez –Q jest zmieniany w zależności od kierunku. Dla –w 2 kolor jest generowany na podstawie kierunku promienia i domyślnego koloru lub zadanego przez –Q (zmiany są gładkie – jest to domyślna opcja). Patrz algorytm koloru tła w 3.5.

-z	Ustawia tryb jednokolorowego raytracingu (3x szybciej). Wynik w odcieniach szarości. Poza tym uruchamia algorytm wykrywający krawędzie trójkątów patrz 3.2.1.6

-Z	Ustawia generowanie oddzielnych obrazów JPEG dla kanałów kolorów: czerwony, zielony, niebieski: powstaną dodatkowe pliki: obraz_r.jpeg, obraz_g.jpeg, obraz_b.jpeg dla kanałów R,G,B.

-F procent	Ustawia algorytm minimalizujący przy generacji drzewa AABB. 0 oznacza algorytm szybki, 100 pełny, wartości (0,100) częściowy, np. –F 8. Domyślnie 100. Podanie wartości 200 lub więcej spowoduje użycie algorytmu “smart minimalize” opisanego w 3.2.1.5 Podanie wartości ujemnych powoduje zastosowanie algorytmu wokselowego, a algorytm minimalizacji wokseli to minus liczba podana, np. Podanie -200 oznacza użycie alorytmu wokselowego oraz algorytm minimalizacji 200 czyli “smart minimalize”, ta opcja jest domyślna. Podanie np -80 oznacza uzycie wokseli, z minimalizacją częściową 80%.  Patrz 3.2.1.

-1 M	ustawia parametry algorytmu sąsiedztwa trójkątów, Liczba M oznacza ilość pikseli generowanych za pomocą jednego drzewa indeksowego, domyslna wartość 1, 0-wyłącza drzewo indeksowe, >1 powoduje znaczne przekłamania wyników, nie zalecane, patrz 3.2.1.6

-H “K N”	Patrz opis w 3.2.1.5. Ustawia wartości N (kiedy dzielić) oraz K (na ile dzielić) algorytmu wokselowego. Domyślne wartości to K=2 i N=7500.

-v liczba	ustawia metodą triangulacji powierzchni NURBS, -v 1 oznacza triangulację równomierną, podanie wartości większej niż 1 oznacza gęstszą triangulację na brzegach a wartości mniejszej niż 1 gęstszą w środku powierzchni, patrz 2.4, np. –v 1.4, domyślnie 1.2, zalecany przedział wartości to <1.,2.>. Uwaga, użycie własnej wartości generuje w wyniku drzewo AABB dla tylko i wyłącznie tej wartości –v, co gorsza ilość trójkątów się nie zmieni więc drzewo dla innej wartości –v też się wczyta, ale boxy tego drzewa będą otaczały inne trójkąty. Można spodziewać się wtedy dziwnych efektów i niezdefiniowanego działania!

-V liczba	skalowanie ilości trójkątów, podanie np. –V 2 oznacza, że chcemy wygenerować w kierunku u i v po 2x tyle trójkątów w stosunku do ilości podanej  w definicji powierzchni NURBS. Zalecany zakres wartości <0.33,3>. Uwaga, AABB drzewa obliczone w preprocesingu tracą swoją ważność – gdyż zmienia się ilość końcowa trójkątów!, np. –V 0.5, domyślnie 1

-j port	na którym porcie uruchomić serwer rays nasłuchujący poleceń. Podanie tej opcji automatycznie uruchamia nowy wątek serwera na zadanym porcie. Porty o numerach mniejszych niż 1024 mogą być niedostępne dla użytkowników nie posiadających uprawnień administratora (nie dotyczy to użytkowników niektórych systemów Windows), np. –j 2500, domyślnie opcja ta jest wyłączona

-B	włącza generowanie sceny binarnej (BIN) na podstawie sceny wejściowej oraz generowanie pliku sceny po przeliczeniu wszystkich przekształceń i triangulacji (plik jest na ogół dużo większy). Plik ten (FDAT) zawiera pełną definicję sceny jako wylistowanie wszystkich trójkątów składowych, nie zawiera powierzchni NURBS ponieważ są w nim zapisane jako trójkąty. Opcja ta jest domyślnie wyłączona

-A	włącza antyaliasing zmniejsza rozdzielczość poziomą i pionową dwukrotnie, ta opcja jest domyślnie wyłączona.

-K	uruchamia opcję zapisu sceny z antyaliasingiem  i bez (2 razy większa rozdzielczość), aby można było kontynuować przerwany raytracing scen z antyaliasingiem włączonym, generuje dodatkowy plik zapisu pełnej sceny z antyaliasingiem np. plik.bmp  plik_aa.bmp, domyślnie ta opcja nie jest aktywna

-2	podwaja rozdzielczość pionową i poziomą, przy jednoczesnym podaniu –A (antyaliasingu) oznacza to, że rozdzielczość pozostanie bez zmian, a włączony zostanie antyaliasing, dodatkowo podanie –K oznacza, że przerwany raytracing z antyaliasingiem , będzie mógł być wznowiony.

-l	wyłącza światło, scena będzie domyślnie maksymalnie oświetlona ze wszystkich stron, domyślnie ta opcja jest nieaktywna

-e	wyłącza teksturowanie, domyślnie teksturowanie jest włączone

-I	wyłącza rozbłyski światła, domyślnie rozbłyski światła są włączone

-u	włącza generowanie cieni, domyślnie cienie nie są generowane

-J 	włącza obsługę JPEG, program wczytując pliki bmp jeżeli nie znajdzie bitmapy to będzie próbował odczytać plik jpeg, np. 1.bmp  1.jpeg, przy zapisie będzie zapisywał zarówno bitmapę jak i plik JPEG, domyślnie opcja ta jest wyłączona

-g	Włącza generowanie JPEG’a w odcieniach szarości, musi być użyta wraz z opcją –J, zapisuje pliki: plik.bmp, plik.jpeg, plik_gs.jpeg

-G 	włącza wyświetlanie w oknie OpenGL w osobnym wątku, patrz: 5.3.3, domyślnie ta opcja jest wyłączona

-f	Uruchamia szybki podgląd sceny w OpenGL, patrz dokładny opis w 5.3.3, ta opcja musi być użyta z opcją –G, domyślnie jest wyłączona

-L	wyłącza reakcje na wszelkie sygnały (SIGUSR1, SIGINT, SIGSEGV), domyślnie program reaguje na te sygnały. SIGUSR1 -> zapisuje obraz na żądanie,  SIGINT zapisuje obraz częściowy (przerwanie), SIGSEGV, zapisuje obraz jak w błędzie krytycznym (panic) i czeka na dołączenie debuggera.

-C	włącza zapisywanie sceny wygenerowanej w preprocesingu do pliku BTREE w formacie tekstowym, np. scene.dat -> scene.btree, domyślnie drzewo nie jest zapisywane do pliku

-E	włącza zapisywanie sceny wygenerowanej w preprocesingu do pliku BTREE w formacie binarnym, np. scene.dat -> scene.btree, domyślnie drzewo nie jest zapisywane do pliku

-c	włącza odczyt preprocessingu sceny z pliku, dla scece.dat próbuje wczytać scene.btree, jak to się nie powiedzie to tworzy drzewo od nowa i jeżeli użyto –C lub –E to ją zapisuje w scene.btree

-N	włącza przetwarzanie powierzchni chropowatych (rozproszenie normalnych), domyślnie rozproszenia nie są liczone

-O	włącza randomizację drzewa w trakcie tworzenia, likwiduje to tendencje niektórych drzew do rozrastania się w określonym kierunku.

-h	wyświetla pomoc programu rays

Powyższa lista nie jest kompletna, aby uzyskać kompletną listę uruchom program rays z opcja –h (pomoc)


6.2.1.2 Opis poleceń w oknie graficznym OpenGL

	Gdy używany OpenGL GUI w trakcie raytracingu (opcja –G) i jednocześnie nie używamy podglądu (opcja –f), to możemy użyć następujących klawiszy.

s	Wysyła sygnał SIGUSR1 do wątku RT, powoduje wygenerowanie obrazu na życzenie
k	Wysyła sygnał SIGINT do wątku RT, powoduje przerwanie działania algorytmu i wygenerowanie obrazu częściowego (działanie jak CTRL+C)
K	Wysyła sygnał SIGKILL bezwarunkowo/natychmiast unicestwiając wszystkie wątki, nic nie jest zapisywane
J/j			manipulacja jakością JPEG’a
A/a			manipulacja wartością oświetlenia tła
M/m			manipulacja wartością minimalnego cienia
V/v			manipulacja wartością maksymalnego cienia
B/b	manipulacja automatycznym zapisem obrazów częściowych zwiększanie/zmniejszanie liczby linii co ile automatyczny zapis
R/r	manipulacja maksymalnym poziomem rekursji
p	przełączanie generacji JPEG’a (włączone/wyłączone)
t			wstrzymaj/wznów raytracing
d			włącz/wyłącz rozpraszanie normalnych
iIuU	zmień współczynniki chropowatości wszystkich trójkątów (i – o mała wartość, u – o dużą wartość)
t	zmień algorytm interpolacji (stary/nowy)
l	zmień algorytm przecięcia (stary bez AABBTree/nowy)
1	włącz/wyłącz światło
2	włącz/wyłącz tekstury
3	włącz/wyłącz rozbłyski
h	wyświetl pomoc

Powyższa lista nie jest kompletna, aby uzyskać kompletną listę uruchom program rays z opcja –h (pomoc)


6.2.1.3 Opis poleceń w oknie podglądu OpenGL

Gdy używamy OpenGL GUI w trybie podglądu (opcja –G –f). Tryb podglądu służy głównie do wygenerowania odpowiedniego przekształcenia sceny, a następnie zapisania go w pliku world_trans.dat za pomocą klawisza w, oto dostępne klawisze w trybie podglądu.


ZzXxCc		Translacje w kierunkach: x,y,z
AaSsDdVv		Skalowanie po: x,y,z, oraz wszystkie osie naraz
123456		Rotacje wg osi x,y,z
lLtTbB		Włącz/wyłącz: światło, teksturowanie, przezroczystość.
UuIiOo		Translacje światła x,y,z
e			Włącz/wyłącz wyświetlanie punktów kontrolnych powierzchni NURBS
7			Włącz/wyłącz podgląd drzewa AABB
890-=			Przemieszczaj się po drzewie AABB
Rr			Przełącz rysowanie: solid/wireframe
n			Odwróć normalne
p			Wyświetl aktualne przekształcenia
w	Zapisz przekształcenia jako ListTransform do pliku world_trans.dat, Uwaga nie jest to kompletny plik DAT, nie należy go wczytywać, należy skopiować z niego przekształcenia do sceny na której działamy za pomocą dowolnego edytora tekstowego, np. Vi.


Powyższa lista nie jest kompletna, aby uzyskać kompletną listę uruchom program rays z opcja –h (pomoc)


6.2.1.4 Opis poleceń serwera rays.


	Aby wydawać polecenia zdalnie należy użyć programu terminal (można też użyć programu telnet, jeżeli użytkownik zna protokół wymiany poleceń z serwerem rays), który służy do wydawania poleceń serwerowi rays, opis terminalu znajduje się w 6.2.2, należy pamiętać o podaniu odpowiedniego portu i adresu IP (lub nazwy) maszyny na której działa proces raytracingu. Można też zdalnie pobrać aktualnie wygenerowany obraz za pomocą programu getbmp opisanego w 6.2.3. Oto dostępne polecenia które można wysyłać do serwera rays:

qserver	zakańcza działanie wątku serwera na zdalnej maszynie, dalsze połączenia będą niemożliwe.
rthlt	wstrzymuje raytracing
rtres	wznawia raytracing
gdist	ustawia globalne rozproszenie normalnej
lidis	wyłącza oświetlenie
tedis	wyłącza teksturowanie
shdis	wyłącza rozbłyski
liena	włącza oświetlenie
teena	włącza teksturowanie
shena	włącza rozbłyski
tmout	zmienia czas odświeżania okna OpenGL
demand	wysyła sygnał SIGUSR1 do wątku RT
intr	wysyła sygnał SIGINT do wątku RT
kill	wysyła sygnał SIGKILL do wątku RT
signal liczba	wysyła zadany sygnał do RT
jqual	manipulacje jakością JPEG’a
tjpeg	włącz/wyłącz generację JPEG’a
tgjpeg	włącz/wyłącz generację JPEG’a w odcieniach szarości
tndist	włącz/wyłącz rozpraszanie normalnych
ambl	manipulacja oświetleniem tła
mins	manipulacja wartością minimalnego cienia
maxs	manipulacja wartością maksymalnego cienia
bkup	manipulacja automatycznym zapisywaniem częściowego obrazu
	(zmniejszanie/zwiększanie ilości linii co ile jest zapis)
recl	zmniejszanie/zwiększanie maksymalnego poziomu rekursji
tgnorm	zmiana algorytmu interpolacji (stary/nowy)
tialg	zmiana algorytmu przecięcia (stary/nowy)
get 	pobiera bitmapę dotychczas wygenerowanego obrazu, nie powinno być używane bezpośrednio, używa program getbmp opisany w 6.2.3
stat	wypisz statystyki działania algorytmu
tex	włącz/wyłącz teksturowanie (tryb podglądu)
light	włącz/wyłącz oświetlenie (tryb podglądu)
tinvn	odwraca normalne (tryb podglądu)
help	wypisuje pomoc

Powyższa lista nie jest kompletna, aby uzyskać kompletną listę uruchom program rays z opcja –h (pomoc)

Przykłady użycia programu RAYS są w rozdziale 7.

6.2.2 Program NURBS

Program ten jest moim projektem z Metod Numerycznych. Generuje on powierzchnię BSPLINE interpolującą zadane punkty. Używa on algorytmu interpolacji globalnej.


	6.2.2.1 Algorytm interpolacji globalnej

ObliczInterpolacje(n,p,s)
Wejście: n – ilość punktów interpolacji, p – wymiar powierzchni, s – wymiar przestrzeni, knot: węzły (knots), t: węzły (nodes), D – tensor punktów interpolacji
(n+1,n+1,s)
Wyjście: powierzchnia BSpline

	OdczytajPunktyInterpolacji

	Pkt_kontrolne = {}
	N = macierz[ n+1]  x [ n+1]
	Policz_N(N)
	NI = OdwróćMacierz(N)

	Policz_Q(Q)
	Policz_P(P)

	Oblicz_Punkty() // wylicza odpowiednio gęstą triangulację powierzchni
			     // korzystając z obliczonych punktów kontrolnych: P

Oblicz_N()
Wejście: macierz N
Wyjście: macierz N z unormowanymi wartościami f-cji bazowych spline

Dla i = 0 do n wykonaj:
	Dla j = 0 do n wykonaj:
		N[i][j] = funkcja_bazowa_bspline(t[i], p, j)
Dla i = 0 do n wykonaj:
	Suma = 0
	Dla j = 0 do n wykonaj suma = suma + N[i][j]
	Dla j = 0 do n wykonaj N[i][j] = N[i][j] / suma

Oblicz_Q()
Wejście: macierz NI, n, s
Wyjście: obliczony tensor Q

Tensor 3-go rzędu: tmp (n+1,n+1,s)
Dla i = 0 do n wykonaj:
	Tmp[i] = NI * D[i]
Dla i = 0 do n wykonaj:
Dla j = 0 do n wykonaj:
Dla k = 0 do s wykonaj: Q[j][i][j] = tmp[i][j][k]

Oblicz_P()
Wejście: tensor Q, macierz NI
Wyjście: tensor punktów kontrolnych P

Dla i = 0 do n wykonaj:	P[i] = NI * Q[i]


6.2.2.2 Opcje programu NURBS

	Ogólna składnia wywołania programu NURBS jest następująca:

NURBS
d|n (l[w|o] file | f[w|o] dim npts def | r[w|o] dim npts | u[w|o] dim npts def1 def2 def3)

Gdzie pierwszy argument d|n (oznaczenie | oznacza, że należy wybrać ‘d’ lub ‘n’). Wybranie ‘d’ spowoduje uruchomienie algorytmu deformacji kostką Beziera, zaś wybranie ‘n’ spowoduje pominięcie deformacji. Zalecaną opcją jest ‘n’ ponieważ deformacje nie dotyczą algorytmu raytracingu i są opcją dostępną w programie NURBS do innych celów (program NURBS jest  projektem zaliczeniowych z przedmiotu Metody Numeryczne i w tej pracy dołączony jest jako użyteczne narzędzie do tworzenia/edycji powierzchni NURBS i tylko w tym celu).

Drugi parametr może mieć jak widać wiele postaci. Ogólnie jest to argument dwuliterowy. Pierwsza litera to: ‘l’ lub ‘f’ lub ‘r’ lub ‘u’, natomiast druga to ‘w’ lub ‘o’. Dozwolone kombinacje są więc następujące: lo, lw, fo, fw, ro, rw, uo, uw.

Pierwsza litera determinuje sposób pozyskania danych:

‘l’:	wczytuje z pliku (wtedy trzeci argument to nazwa pliku). Pliki można zapisywać w trakcie działania programu z menu kontekstowego
‘f’:	nakazuje użyć funkcji jako źródła danych, wtedy następne trzy parametry to kolejno: wymiar powierzchni spline, ilość punktów i definicja funkcji. Wymiar powinien być co najmniej 2 aby funkcja była gładka, ilość punktów musi być co najmniej wymiar + 1. Jeżeli podamy wymiar np.: 3 a ilość punktów 5, to funkcja zostanie spróbkowana w 5 x 5 = 25 punktach (na siatce). Punkty te staną się punktami interpolacji. Ostatni parametr to definicja funkcji. Program jest wyposażony w parser więc można ją podać jako argument, np. ‘sin(x+y)-cos(x-y)’. Czyli cała linia poleceń byłaby np. taka: n fw 3 5 ‘sin(x+y)’
‘r’:	nakazuje użyć losowych danych, należy następnie podać wymiar i ilość punktów. Ilość punktów powinna być >= wymiar + 1, wymiar powinien być >= 2 aby powierzchnia była gładka. Przykładowa linia poleceń: n rw 3 10
‘u’: Nakazuje użyć 3 funkcji dla każdej współrzędnej obrazu, czyli tak jak w ‘f’ ale należy podać 3 definicje funkcji. Podanie ‘x’ i ‘y’ jako dwóch ostatnich funkcji daje rezultat identyczny jak użycie ‘f’. Przykładowa linia poleceń: uw 3 10 ‘sin(x*y) ‘x^2’ ‘y^2’.

Druga litera decyduje o sposobie wyświetlania wyników.
‘w’:	Wyświetlanie poszczególnych trójkątów (zalecane)
‘o’:	Wyświetlanie zoptymalizowane całych powierzchni (przez OpenGL) – ładnie wygląda ale nie zalecane, ponieważ chcemy mieć kontrolę nad triangulacją.

Kilka przykładów uruchomienia:

NURBS n lw bspline_surface1
NURBS n rw 3 10
NURBS n rw 7 8
NURBS n fw 2 15 ’sin(x+y)’
NURBS n fw 4 8 ’x*y’
NURBS n uw 3 10 ’sin(x*y) ’x^2’ ’y^2’
6.2.2.3 Używanie programu NURBS

Program można kontrolować na dwa sposoby: menu kontekstowe oraz klawiatura.

Opiszę teraz opcje dostępne przy zastosowaniu menu kontekstowego:

‘write bspline surface’	zapisuje punkty kontrolne powierzchni w formacie wewnętrznym programu oraz w formacie NURBS
	pliki bspline_surfaceN oraz surface.NURBS
‘write Rayslib NURBS surface’	zapisuje aktualną triangulację w raylab.dat oraz powierzchnię w formacie NURBS w NURBSsurfaceN.dat
‘toggle D/P/Q points’	włącza i wyłącza wyświetlanie punktów interpolacji/pośrednich/kontrolnych
‘more/less u/v lines’			zwiększa/zmniejsza gęstość triangulacji w kierunkach u/v
‘lower/higher degree’	zmniejsza/zwiększa stopień powierzchni, minimalny: 1 maksymalny: ilość punktów – 1.

Opis klawiatury:

1,2,3			włącza i wyłącza wyświetlanie punktów								interpolacji/pośrednich/kontrolnych
ESC,q			zakańcza działanie aplikacji
`	zapisuje aktualną triangulację w raylab.dat oraz powierzchnię w formacie NURBS w surface.NURBS
r	zapisuje punkty kontrolne powierzchni w formacie wewnętrznym programu oraz w formacie NURBS, pliki: bspline_surfaceN oraz NURBSsurfaceN.dat
h	wyświetla pomoc
wsadex	przemieszcza punkt kontrolny (aktualnie zaznaczony)
ikjlom	obraca cały obiekt dookoła różnych osi
p	włącza/wyłącza przemieszczanie aktualnie zaznaczonego punktu
nbvc	przemieszczanie się do następnego/poprzedniego/prawego/lewego punktu kontrolnego
4	włącza/wyłącza automatyczne obroty
+-=	skalowanie
90	zmniejsz/zwiększ stopień powierzchni
5678	zmniejsz/zwiększ gęstość triangulacji u/v
SPACJA	powróć do domyślnych ustawień


	6.2.3 Program NURBS2DAT

	Program służy do konwersji plików NURBS do formatu DAT. Dawna wersja RAYS nie obsługiwała plików NURBS bezpośrednio więc taka konwersja była potrzeba, teraz program RAYS potrafi już bezpośrednio odczytać definicję NURBS w pliku DAT i ztriangulować tę powierzchnię „w locie”.

	Parametry programu, kolejno:
	-plik wejściowy np. powierzchnia.NURBS
	-plik wyjściowy np. scena.DAT
	-index: (od jakiego numeru zaczynać numerację trójkątów) np. 0
	-czy zapisać nagłówek: 0 lub jeden

Przykładowe użycie programu: NURBS2DAT input.NURBS output.DAT 0 1

	6.2.4 Program IGES2DAT

	Program służy do odczytu powierzchni NURBS zapisanych w pliku IGS, a następnie zapisaniu ich w pliku DAT. Program odczytuje wszystkie encje o numerze 128 (powierzchnie NURBS). Program jest bazowany na aplikacji renderującej NURBS’y z plików IGES – projekt z Laboratoriów Systemów CAD/CAM – xcam.

Program obsługuje następujące opcje z linii poleceń:
-h			wyświetla pomoc
-T ‘x y z’		translacja o x,y,z
-R ‘rx ry rz’		rotacja o rx wokoło 0X, ry wokoło 0Y i rz wokoło 0Z
-X ‘x y z’		skalowanie x,y,z
-T ‘r g b’		przezroczystość: składowe koloru RGB
-C ’r g b’		kolor: składowe koloru RGB
-S ’r g b’		odbicie światła: składowe koloru RGB
-d			włącza tryb DEBUG
-N percent		ustawia chropowatość materiału (normal distorber)
-I liczba		ustawia współczynnik rozbłysku dla wszystkich 					kolorów
-m tid			ustawia index textury na tid
-t liczba		gęstość triangulacji liczba * liczba * 2, zalecane: 10-64
-i plik.iges		plik wejściowy
-o plik.dat		plik wyjściowy

Przykładowe użycie programu:
	IGES2DAT –i scena.igs –o scena.dat –X ‘5 5 5’ –T ‘1 1 1’ –S ‘1 1 1’


	6.2.5 Program ULI2DAT

	Program konwertuje pliki ULI/TRI (nasz format na VR’ach) na pliki DAT. Należy podać następujące parametry z linii poleceń:
index	początkowy index numeracji. Jeżeli np. w pliku ULI będzie 1000 trójkątów i podamy index 100 to w pliku DAT zostaną zapisane od 100 do 1100
tr tg tb	Współczynniki przezroczystości dla kolejnych kolorów RGB np        0.2 0.4 0.35
sr sg sb	Współczynniki odbicia światła dla kolejnych kolorów RGB np.
	1	0.6 0.85
cr cg cb	Kolor materiału (RGB) np. 1 0 0 daje czerwony, 0 1 0 zielony, np.
		0.5 0.5 0.1
tfr tfg tfb	Współczynniki załamania światła dla kolejnych kolorów np.
		1.22 1.25 1.27
sf		Współczynnik rozbłysku dla wszystkich kolorów, np. 64.0
tid		index tekstury
fac	1 lub 2, czy materiał jest jednostronny czy dwustronny. Na dwustronnych materiałach rozbłysk i oświetlenie jest po obu stronach
rand	0 lub 1, jeżeli włączone to spowoduje użycie losowych wartości
hdr	0 lub 1, 1-zapisuje nagłówek pliku, 0-nie zapisuje

Przykładowe użycie programu:
ULI2DAT 0 0.25 0.25 0.25 0.5 0.5 0.5 0.25 0.25 0.25 1.22 1.25 1.27 64.0 2 2 0 1 < input.ULI > output.DAT


	6.2.6 Program MD22DAT

	Program konwertuje pliku MD2 (modele z gry Quake II) do formatu DAT. Należy podać następujące parametry z linii poleceń:

tfr tfg tfb	Współczynniki załamania światła dla kolejnych kolorów np.
			1.22 1.25 1.27
tr tg tb	Współczynniki przezroczystości dla kolejnych kolorów RGB np        0.2 0.4 0.35
sr sg sb	Współczynniki odbicia światła dla kolejnych kolorów RGB np.
	1	0.6 0.85
cr cg cb	Kolor materiału (RGB) np. 1 0 0 daje czerwony, 0 1 0 zielony, np.
		0.5 0.5 0.1
sf		Współczynnik rozbłysku dla wszystkich kolorów, np. 64.0
nd		Współczynnik chropowatości (normal distorber)
fac	1 lub 2, czy materiał jest jednostronny czy dwustronny. Na dwustronnych materiałach rozbłysk i oświetlenie jest po obu stronach
tid	index tekstury
index	początkowy index numeracji. Jeżeli np. w pliku MD2 będzie 1000 trójkątów i podamy index 100 to w pliku DAT zostaną zapisane od 100 do 1100
hdr		0 lub 1, 1-zapisuje nagłówek pliku, 0-nie zapisuje

Obsługa klawiatury:

ESC/q	zakończenia działania programu
SPACJA	wykonaj konwersję modelu (aktualnie wyświetlanej ramki animacji)
h		pomoc
d		włącz/wyłącz generowanie nagłówka pliku DAT
-=		zmniejsz/zwiększ szybkość animacji
xXyYzZ	rotacje wokół x/y/z
12		skalowanie

Przykładowe użycie programu:
MD22DAT 1.22 1.25 1.27 1 1 1 2 2 2 1 1 1 64 2.5 2 1 100 1 < ogre.MD2 > ogre.DAT

	6.2.7 Program 3DS2TRI

	Program służy do konwersji obiektów zapisanych w formacie 3DS do formatu TRI/ULI, który może być następnie skonwertowany do formatu DAT za pomocą programu ULI2DAT. W trakcie konwersji tracone są dane o teksturach. Ponieważ istnieje już wersja konwertująca bezpośrednio z 3DS do DAT (i zachowująca tekstury), więc program ten nie jest już na ogół używany. Opcje tego programu są podobne jak w programie 3DS2DAT.
	6.2.8 Program 3DS2DAT

	Program konwertuje pliki 3DS do formatu DAT. Normalne (których nie ma w obiekcie 3DS są interpolowane dla każdego wierzchołka, dzięki czemu powstaje obiekt składający się z trójkątów GPT. Brak interpolacji powodował, zę obiekty wyglądały nienaturalnie, odpowiednie obrazki są załączone na koncu pracy. Można zapisać przekształcenia obiektu, nadać mu właściwości materiałowe itp. Oto dostępne opcje z linii poleceń, każdą opcję podaje się w formacie ‘opcja=wartość’:

‘cr=liczba’			wartość koloru czerwonego
‘cg=liczba’			wartość koloru zielonego
‘cb=liczba’			wartość koloru niebieskiego
‘tr=liczba’			przezroczystość koloru czerwonego
‘tg=liczba’			przezroczystość koloru zielonego
‘tb=liczba’			przezroczystość koloru niebieskiego
‘sr=liczba’			odbijalność koloru czerwonego
‘sg=liczba’			odbijalność koloru zielonego
‘sb=liczba’			odbijalność koloru niebieskiego
‘tfr=liczba’			współczynnik załamania koloru czerwonego
‘tfg=liczba’			współczynnik załamania koloru zielonego
‘tfb=liczba’			współczynnik załamania koloru niebieskiego
‘faces=liczba’		1 lub 2 (patrz opis fac w 6.2.6)
‘tx=liczba’			translacja w kierunku x
‘ty=liczba’			translacja w kierunku y
‘tz=liczba’			translacja w kierunku z
‘rx=liczba’			rotacja dookoła osi 0X
‘ry=liczba’			rotacja dookoła osi 0Y
‘rz=liczba’			rotacja dookoła osi 0Z
‘scale=liczba’		skalowanie
‘notex’			pomiń tekstury

	Obsługa klawiatury:

q				zakończ działanie aplikacji
n				włącz/wyłącz interpolowanie normalnych w wierzchołkach
wsadex			rotacje
ikjlom			translacje
yu				skalowanie
12				zapisz trójkąty binarnie/tekstowo
3				zapisz scenę jako plik DAT: model.3ds  model.3ds.DAT

Przykładowe użycie programu:
3DA2DAT model.3ds ‘cr=1’ ‘cg=1’ ‘cb=1’ ‘sclae=10’


	6.2.9 Program 60FACES

	Program generuje dwunastościan foremny składający się z 60 trójkątów (po 5 trójkątów na każdą ścianę). Należy podać następujące parametry z linii poleceń:

	Index tr tg tb sr sg sb cr cg cb tfr tfg tfb zraise sf tid fac lt

Większość z tych parametrów została już opisana powyżej, opiszę tylko te które nie były dotychczas opisywane:

	zraise: o ile podnieść wierzchołek w środku każdego pięciokąta w stosunku do płaszczyzny ściany. Podanie 0 wygeneruje 12-ścian foremny, inne wartości spowodują wygenerowanie 60-ściana. Podanie wartości 100 spowoduje wybranie takiej wartości dla której wszystkie trójkąty będą równoboczne

	lt: podanie 1 wygeneruje tylko listę przekształceń, podanie 0 wygeneruje tylko trójkąty, podanie 2 spowoduje wygenerowanie obu składników

	Efektem działania programu jest fragment pliku DAT który może być użyty do wstawienia 12/60-ścianów do dowolnej sceny (generowane może być też ListTransform dla danego obiektu dla lt=1 lub 2)

	Przykładowe użycie programu:
		60FACES 1 1 1 0 0 0 1 1 1 1.1 1.2 1.3 0 48.0 2 1 2 > fragment.DAT


	6.2.10 Program TORRUSGEN

	Program generuje torrus na podstawie zadanych parametrów. Powstaje fragment pliku DAT który może być użyty w innych scenach. Należy podać następujące parametry z linii poleceń:

R r N n index cr cg cb sr sg sb tr tg tb tid faces sf tfr tfg tfb

Większość z tych parametrów została już opisana powyżej, opiszę tylko te które nie były dotychczas opisywane:


	R r: promienie torrusa: duży i mały
	N n: ilość podziałów u i v torrusa, ilość trójkątów = 2*N*n


Przykładowe użycie programu:
	TORRUSGEN 100 25 24 24 0 1 0 0   0 1 0    0 0 1 10 1 128.0 1.25 1.26 1.27 > t.DAT


	6.2.11 Program CONE

	Program generuje stożek na podstawie zadanych parametrów. Powstaje fragment pliku DAT który może być użyty w innych scenach. Należy podać następujące parametry z linii poleceń:

r h n index cr cg cb sr sg sb tr tg tb tid faces sf tfr tfg tfb

Większość z tych parametrów została już opisana powyżej, opiszę tylko te, które nie były dotychczas opisywane:

	r: promień stożka
	h: wysokość stożka
	n: ilość podziałów podstawy stożka: Ilość trójkątów = 2*n


Przykładowe użycie programu:

CONE 10 30 120  0 1 0 0   0 1 0    0 0 1 10 1 128.0 1.25 1.26 1.27 > c.DAT

	6.2.12 Program CUBE

	Program generuje sześcian (12 trójkątów) na podstawie zadanych parametrów. Powstaje fragment pliku DAT który może być użyty w innych scenach. Należy podać następujące parametry z linii poleceń:

index tr tg tb sr sg sb cr cg cb tid tfr tfg tfb sf sx sy sz rx ry rz tz ty tz invN fac lt

Większość z tych parametrów została już opisana powyżej, opiszę tylko te które nie były dotychczas opisywane:

	sx sy sz: skalowanie sześcianu
	rx ry rz: rotacje sześcianu
	tx ty tz: translacje sześcianu
	invN: 1-odwraca normalne do wewnątrz, 0-normalne bez zmian
	lt: 0,1,2 – działanie jak w 6.2.9.


Przykładowe użycie programu:

CUBE 0 1 1 1 1 1 1 1 1 1 2 1.2 1.3 1.4 64.0 4 2 3 30 20 10 –10 –20 50 0 1 2 > cube.DAT

	6.2.13 Program CYLINDER

	Program generuje walec na podstawie zadanych parametrów. Powstaje fragment pliku DAT który może być użyty w innych scenach. Należy podać następujące parametry z linii poleceń:

r h n index cr cg cb sr sg sb tr tg tb tid faces sf tfr tfg tfb

Większość z tych parametrów została już opisana powyżej, opiszę tylko te które nie były dotychczas opisywane:

	r: promień walca
	h: wysokość walca
	n: ilość podziałów podstawy walca: Ilość trójkątów = 4*n


Przykładowe użycie programu:

CYLINDER 10 30 120  0 1 0 0   0 1 0    0 0 1 10 1 128.0 1.25 1.26 1.27 > cyl.DAT

	6.2.14 Program BALL

	Program generuje kulę jednostkową na podstawie zadanych parametrów. Powstała kula składa się z dwóch niezależnych półkul, dla których są generowane osobne listy przekształceń. Powstaje fragment pliku DAT który może być użyty w innych scenach. Należy podać następujące parametry z linii poleceń:

index tr tg tb sr sg sb cr cg cb tid tfr tfg tfb sf sx sy sz rx ry rz tz ty tz invN fac lt pa pb

Większość z tych parametrów została już opisana powyżej, opiszę tylko te które nie były dotychczas opisywane:


	pa: ilość podziałów parametru alfa
	pb: ilość podziałów parametru beta, ilość trójkątów 4*pa*pb


Przykładowe użycie programu:

BALL 0 1 1 1 1 1 1 1 1 1 2 1.2 1.3 1.4 64.0 4 2 3 30 20 10 –10 –20 50 0 1 2 16 16 > b.DAT


	6.2.15 Program RANDNURB/RANDNURBFULL
Oba programy przyjmują tylko jeden argument z linii poleceń: ilość powierzchni NURBS do wygenerowania. Generują one losowe powierzchnie NURBS (program RANDNURB generuje bardziej płaskie i stabilne powierzchnie, a RANDNURBFULL generuje bardzo losowe powierzchnie). Oba programy zapisują na stdout, wynik pliku można zapisać jako plik NURBS używając przekierowania.

Przykładowe użycie programu:

RANDNURB 10 > 10surfaces.NURBS
RANDNURBFULL 1000 > 1000surfaces.NURBS

	6.2.16 Program RTRIANGLE

	Program przyjmuje jeden argument z linii poleceń: ilość trójkątów do wygenerowania. Wynik jest wypisywany na stdout, można go zapisać w pliku DAT jako samodzielną scenę.

	Przykładowe użycie programu:

	RTRIANGLE 100 > scena.dat

	6.2.17 Program TABLE

	Program generuje stół używając do tego polecenia CUBE, wszystkie parametry podaje się z linii poleceń w następującej kolejności:

	Index tid1 tid2 tid3 size cubecmd

	Index:		indeks od którego zaczną się trójkąty stołu
	Tid1,2,3:		tekstury kolejnych elementów stołu: noga, blat, szkło
				sugerowane: 8 4 9
	Size:			wielkość stołu 9skala)
	Cubecmd:		polecenie które należy wywołać aby wygenerować sześcian:
				Np.: CUBE, cube, cyg_cube.exe itp.

	Jak widać ten program jest właściwie skryptem generującym odpowiedni zestaw prostopadłościanów tworzących stół.

	Przykładowe użycie programu:

	TABLE 0 8 4 9 10.5 ‘./cube’


	6.2.18 Program TERMINAL


	Program służy do nawiązywania połączeń z serwerem RAYS oraz do wydawania poleceń serwerowi. Opcje programu są następujące:


-i  ‘adresIP’		podaj adres IP serwera, domyślny ‘127.0.0.1’
-s  serwer		podaj nazwę serwera, obliczy IP używając systemowego DNS
-p port	podaj port na którym nasłuchuje serwer, domyślny 2500. Jest to numer portu, który został podany za pomocą opcji –j w programie RAYS (patrz 6.2.1.1)
-c polecenie	polecenie wysyłane do serwera RAYS, lista dostępnych poleceń znajduje się w 6.2.1.4.

Przykładowe użycie programu:

TERMINAL –i ‘213.100.91.147’ –p 1999 –c ‘stat’


	6.2.19 Program GETBMP

	Program służy do pobrania aktualnie obliczonej bitmapy z serwera RAYS. Opcje programu są następujące:


-i  ‘adresIP’		podaj adres IP serwera, domyślny ‘127.0.0.1’
-s  serwer		podaj nazwę serwera, obliczy IP używając systemowego DNS
-p port	podaj port na którym nasłuchuje serwer, domyślny 2500. Jest to numer portu, który został podany za pomocą opcji –j w programie RAYS (patrz 6.2.1.1)
-o plik.bmp	nazwa pliku wynikowego, domyślnie „output.bmp”


Przykładowe użycie programu:

GETBMP –s starlight64 –p 1250 –o map.bmp
	6.2.20 Program TEX

	Program  generuje przykładową teksturę w formacie BMP. Tekstura jest samopowtarzalna. Należy podać następujące parametry z linii poleceń:

	Texname: nazwa tekstury wynikowej, np.: tex.bmp
Współczynnik potęgowy: liczba kontrolująca algorytm wyliczający kolory, zalecane: 0.2-16, im wyższa liczba tym „ostrzejsze” zmiany koloru

Przykładowe użycie programu:

	TEX tex1.bmp 8.0

	6.2.21 Program WRAPPER

	Program uruchamiający RAYS z odpowiednimi parametrami. Jest to interaktywna nakładka tworząca skrypt wywołujący RAYS, pyta się o wszystkie opcje użytkownika i zapisuje ostateczny wiersz poleceń w pliku: scripts/rays_cmdN.sh. Program nie oczekuje żadnych opcji z linii poleceń, należy odpowiadać na pytania zadawane przez program:

	Użycie programu:

	WRAPPER

	6.2.22 Program BTREECONV

	Jest to program do manipulacji AABB-drzewami. Obsługiwane opcje to:
	-odczyt drzewa AABB w formacie binarnym lub tekstowym
	-zapis drzewa AABB w formacie binarnym lub tekstowym
	-dokonywanie transformacji drzew (skalowania, translacje i niektóre rotacje)
	-łączenie dwóch drzew w jedno.

	Transformacje wykonywane są w następującej kolejności: skalowanie, rotacja, translacja. Ponieważ drzewa AABB z definicji muszą być wyrównane do osi układu współrzędnych to bez regeneracji drzewa można wykonać tylko następujące obroty: o 90, 180 i 270 stopni dookoła osi OX, OY lub OZ. Tylko takie obroty są obsługiwane, ponieważ program nie służy do regeneracji drzewa (generacją drzewa zajmuje się RAYS) . Program służy do szybkich modyfikacji AABB drzew nie wymagających ich tworzenia od nowa (co potrafi trwać godzinami przy dużych obiektach, a przekształcenie zajmuje najwyżej kilka sekund).

	Przy łączeniu drzew należy pamiętać o następujących rzeczach. Po pierwsze najlepiej łączyć drzewa obiektów rozłącznych – wtedy otrzymujemy w pełni optymalne drzewo, bardzo małym kosztem. np. generujemy drzewo dla obiektu mającego 40K trójkątów (około kilku godzin), następnie przekształcamy to drzewo translacją o 100 w kierunku np. X (około 2-3 sekund), kopiujemy drzewo systemowym poleceniem i przekształcamy o np -200. Mamy dwa drzewa, teraz łączymy je (około 2-3 sekund) i mamy drzewo dwóch obiektów (najlepiej aby były one rozłączne, tj ich wielkości w kierunku X nie były większe niż 100). W ten sposób mamy drzewo 80K w kilka sekund. Bez programu do konwersji drzew musielibyśmy stworzyć scenę z 80K trójkątów i dla niej wygenerować drzewo w całości (nie korzystając z rozłączności), potrwałoby to około doby!, ponieważ złożoność generatora drzew to O(n^3). Inny przykład: mamy obiekt mający 10K trójkątów, czas generacji około 1 minuty. Kopiujemy go 16 razy, mamy obiekt 160K czas generacji drzewa około 1-2 tygodnie!. A używając kopiowania drzewa obiektu to: 1min wygenerowanie jednego drzewa, a następnie kopiowanie 15 razy, powiedzmy około 15-30 min, zysk jest duży.

	Jeżeli drzewa łączone przecinają się to drzewo wynikowe będzie poprawne, ale nie optymalne. Jednak strata optymalności jest niewielka. Dla drzewa oddzielnych obiektów powstanie po prostu korzeń decyzyjny czy pójść do jednego obiektu czy do drugiego i dalej odpowiednio drzewo 1 lub 2 będzie przetwarzane. Dla drzew przecinających okaże się, że dla pewnych promieni będziemy sprawdzać oba drzewa, ale strata będzie niewielka, poza tym gdybyśmy wygenerowali całe drzewo bez łączenia to też byśmy musieli przeglądać elementy jednego i drugiego obiektu. Trudno oszacować stratę optymalności ale dla testów na łączeniu obiektów 4K i 4K nie zauważyłem istotnej straty wydajności: wysokości obu drzew (łączonego i generowanego) były identyczne, struktura się różniła, średnia ilość przecięć była podobna Dla rozłącznych obiektów otrzymałem: Ilość przecięć drzewa/trójkątów dla drzewa łączonego programem btreeconv:
	3541064/958266 co daje odp procenty 0.227%/0.061%
	A dla wygenerowanego dla obu obiektów
	3635896/958012 co daje odp procenty 0.233%/0.061%
	procenty oznaczają ile było procent przecięć drzewa/trójkątów w stosunku do liczenia bez lokalizacji (brutal force). Jak widać dla rozłącznych obiektów przecięć z drzewem było nawet mniej niż w generowanym drzewie (o około 89000) za to ilość przecięć trójkątów była większa ale bardzo nieznacznie: o 254 przecięcia.
	Sprawa ma się gorzej dla przecinających się obiektów (utaj kule o promieniu 100 ze środkami w -10,0,0 i 10,0,0). Oto wyniki dla łączenia drzew:
	3181536/754462 procenty: 0.245%/0.058%
	A dla generowanego dla obu obiektów naraz (i bez łączenia):
	2610723/777926 procenty: 0.201%/0.060%
	Jak widać drzewo generowane dla obu obiektów jest rzadziej przecinane (jest bardziej dopasowane/optymalne), zaś drzewo łączone jest przecinane dużo częściej (około 20%)  jednak więcej przecięć drzewa nie generuje przecięcia z trójkątem.
	Dla przypadku połączenia drzew dla identycznych obiektów mamy
	3422327/930272 procenty: 0.271%/0.074%
	A dla wygenerowania dla dwóch identycznych obiektów
	2152985/737084 procenty: 0.171%/0.058%
	Jak widać jest to najgorszy możliwy przypadek, dodatkowy nakład na przecinanie drzewa wynosi ponad 50%, zaś przecięć trójkątów jest o około 25% więcej. Nie jest to zły wynik.

	Należy także pamiętać o poprawnej kolejności łączenia drzew obiektów (zachowując jak najbardziej drzewiasta strukturę). Np. Można tak połączyć obiekt: a=b+c, d=a+e, f=d+g, h=f+i, j=h+k, l=j+m, n=l+o. otrzymamy w ten sposób strukturę o wysokości 7, a można te obiekty połączyć tak: a=b+c, b=d+e, c=f+g, d=h+i, e=j+k, f=l+m, g=n+o i otrzymamy zrównoważone drzewo o wysokości 3. W ten sposób otrzymamy bardziej optymalne drzewo bo niższe, obliczenia na tym drzewie będą przebiegały szybciej.

Program można uruchomić w dwóch trybach. Tryb transformacji drzewa i tryb łączenia drzew. Aby uruchomić w trybie transformacji drzewa należy podać następujące parametry z linii poleceń:

t|b input.btree output.btree tx ty tz sx sy sz rot

t|b oznacza, że należy podać format pliku wyjściowego: t-tekstowy, b-binarny, formatu pliku wejściowego nie należy podawać, zostaje on automatycznie wykryty.
Input.btree i output.btree to kolejno wejściowy plik drzewa i wyjściowy plik drzewa
koleje 6 parametrów to 3 parametry translacji i 3 parametry skalowania
rot to rotacja drzewa, dozwolone wartości to:
	r0 (brak rotacji)
	x90, x180, x270 (rotacja dookoła OX o 90,180 lub 270 stopni)
	y90, y180, y270 (rotacja dookoła OY o 90,180 lub 270 stopni)
	z90, z180, z270 (rotacja dookoła OZ o 90,180 lub 270 stopni)

Aby uruchomić program w trybie łączenia drzew należy podać następujące parametry z linii poleceń

	mt|mb input1.btree input2.btree output.btree

	mt|mb oznacza, że należy podać format pliku wyjściowego: mt-tekstowy, mb-binarny, formatu plików wejściowych nie należy podawać, zostaną one automatycznie wykryte.
	Kolejne 3 argumenty to nazwy plików, odpowiednio dwóch plików wejściowych i plik wyjściowy będący ich sumą (kolejność jest istotna, taka sama musi być w pliku DAT – tj. Ten drugi obiekt musi być także drugim w pliku DAT). Ostatecznie:
output.btree = input1.btree + input2.btree

Przykładowe użycie programu:
Zamiana z tekstowego na binarny:
	BTREECONV b in.btree out.btree 0 0 0 1 1 1 r0
Transformacje:
	BTREECONV t in.btree out.btree 0 100 0 1 1 1 r0
	BTREECONV b in.btree out.btree 0 0 0 2 2 2 r0
	BTREECONV t in.btree out.btree 0 0 0 1 1 1 z180
Łączenie drzew:
	BTREECONV mb in1.btree in2.btree out.btree


7. Przykładowe sceny, sposób renderowania i metody kompilacji

7.1 Przykładowe sceny

	Przykładowe sceny znajdują się w katalogu DAT/, są to wszystkie sceny stworzone w trakcie tworzenia aplikacji RAYS. Gotowe przykłady demonstracyjne znajdują się w katalogu EXAMPLES/ - są to uruchamialne skrypty zawierające wszystkie opcje wywołania programu RAYS. Kolejne przykłady to exN.sh, gdzie n = {1,2,3....}. Każdy przykład wyświetla najpierw opis, a potem dokonuje obliczenia sceny. Przykłady demonstrują różne opcje programów,przed wykonaniem jakiegokolwiek polecenia wyświetlają to polecenie.


7.2 Kilka przykładów renderowania

RAYS -i nazwapliku.DAT
	– najprostsze renderowanie, brak obsługi JPEG, braku GUI, brak zapisu AABB drzewa, plik wynikowy: screen.BMP
RAYS -i nazwapliku.DAT -J
	⁃	dodanie obsługi tekstur w formacie JPEG, oraz zapisu wyniku w formacie JPEG
RAYS -i nazwapliku.DAT -J -G
	⁃	dodanie wyświetlania w OpenGL
RAYS -i nazwapliku.DAT -J -G -c
	⁃	dodanie obsługi odczytu AABB-drzewa z pliku
RAYS -i nazwapliku.DAT -J -G -C
	⁃	zapisanie wygenerowanego AABB-drzewa
RAYS -i nazwapliku.DAT -J -G -C -c
	⁃	obsługa zapisu/odczytu AABB-drzewa
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10
	- zmiana pliku wynikowego, oraz maksymalna rekursja: 10
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -N
	- obsługa chropowatych powierzchni
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -T /usr/textures
	⁃	zmiana katalogu w którym szukane są tekstury
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -R partial.BMP
	⁃	kontynuacja renderingu pliku partial.BMP
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -x 100 -y 100
	⁃	zmiana rozdzielczości obrazu
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -m 0.2 -M 0.4
	⁃	ustawienie minimalnego cienia na 0.2 a maksymalnego na 0.4
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -a 0.5
	⁃	ustawienie światła tła na 0.5.
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -Z
	⁃	oddzielne obrazy dla kanałów RGB
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -Q ff0000
	⁃	ustawienie koloru tła na FF0000 (czerwony)
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -W 00ff00
	⁃	ustawienie koloru światła na 00FF00 (zielony)
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -z
	⁃	rendering w odcieniach szarości
RAYS -i nazwapliku.DAT -J -G -C -o wynik.BMP -r 10 -F 8
	⁃	ustawia minimalizację na 8% przeszukiwania drzewa
RAYS -i nazwapliku.DAT -J -G -C -o wynik.BMP -r 10 -F 0
	⁃	ustawia minimalizację na najszybszy (i najgorszy) algorytm o(n^2)
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -A
	⁃	włącza antyaliasing, rozdzielczość zmniejsza się dwukrotnie
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -A -2
	⁃	podwaja rozdzielczość w trybie antyaliasingu, więc ostateczna rozdzielczość jest bez zmian
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP MP -r 10 -A -K -2
	⁃	umożliwia kontynuację raytracingu z antyaliasingiem (zapisuje dodatkowo plik o podwojonej rozdzielczości bez antyaliasingu)
RAYS -i nazwapliku.DAT -G -C -c -o wynik.BMP -r 10 -e
	⁃	wyłącza teksturing, wszelkie tekstury i odwołania do nich są pomijane
RAYS -i nazwapliku.DAT -G -C -c -o wynik.BMP -r 10 -j 20000
	- uruchamia serwer rays na porcie 20000
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -u
	⁃	wyłącza obliczanie cieni
RAYS -i nazwapliku.DAT -J -G -C -c -o wynik.BMP -r 10 -I
	⁃	wyłącza rozbłyski
RAYS -i nazwapliku.DAT -J -G -f
	⁃	uruchamia tryb szybkiego podglądu sceny
RAYS -i nazwapliku.DAT -J -G -H “3 5000” -F -200
	⁃	uruchamia algorytm wokselowy z K=3, N=5000 i pełną minimalizacją wokseli

	7.3 Sposoby kompilacji

	Kompilacja w środowisku UNIX'owym odbywa się przez wydanie polecenia make lub gmake (GNU-Make). Program jest napisany w języku C, używa bibliotek: OpenGL, GLUT, JPEGlib. Można także kompilować wersję pozbawione pewnych funkcjonalności, należy dodać odpowiednie flagi kompilacji (patrz plik Makefile). Dostępne flagi to:

	-DNOINET:		nie wkompilowuje obsługi przez internet
	-DNOJPEG:		nie wkompilowuje obsługi JPEG
	-DNOGL:		nie wkompilowuje obsługi OpenGL GLUT
	-DNOSIGNALS:	nie wkompilowuje obsługi sygnałow

	W przypadku braku programu make (głównie systemy Microsoftu), można skorzystać ze skryptów kompilacyjnych: bat/compile_rays.bat i bat/compile_win32.bat
	Zalecany kompilator to GNU GCC 3.X, 2.X powinien działać. Kompilator Microsoftu MSVC NIE jest obsługiwany (i nie będzie). Jedyna możliwość kompilacji w środowisku Windows to użycie pakietu Cygwin. Program testowano na następujących maszynach:

	FreeBSD UNIX 6.0-BETA1 amd64:			platforma 			macierzysta
	FreeBSD UNIX 4.X,5.X, 6.X i386:			działa
	NetBSD UNIX 3.0 amd64:					działa
	NetBSD UNIX 3.0 i386:					działa
	Slackware Linux 10.0 i386:				działa
	SUN Solaris 10 amd64 + GNUtools:			działa
	SUN Solaris 9 i386 + GNUtools:				działa
	SUN Solaris 9 SPARC 					nie działa *
	MS Windows Server 2003 x86 + Cygwin:		działa
	MS Windows XP Home x86 + Cygwin:			działa
	MS Windows XP x64 edition + Cygwin (32bit):	działa
	MS DOS + bcc:						nie działa **

* format bajtów BIGENDIAN, problemy z zapisem/odczytem danych
	** segment kodu/danych przekracza 64K

		7.4 Struktura katalogów

		3ds2dat/, 3ds2tri/:	-programy do konwersji z formatu 3DS, modele (katalog 						3ds2tri/3ds/)
		bat/			- skrypty do kompilacji na systemie Windows
		bmp/			- bitmapy, jpeg'i wygenerowane przez program
		dat/			- pliki wejściowe (sceny)
		iges/			- pliki IGS zawierające przykładowe pow. NURBS
		libjpeg/		- biblioteka jpeg (źródła)
		md2/			- modele z gry Quake II
		pcx/			- tekstury modeli z gry Quake II
		textures/		- domyślny katalog gdzie są poszukiwane tekstury
		win/			- zestaw bibliotek/narzędzi dla systemu Windows

		7.5 Przykładowe obrazy wygenerowane za pomocą RAYS

		Efekt soczewki i rozproszenia światła

Skomplikowana powierzchnia NURBS z teksturą












Pryzmat















Kilka renderingów Sphinx'a składającego się z 58K trójkątów..








Obiekt wczytany z pliku 3DS...





Soczewka jako powierzchnia NURBS...


















Inne obiekty 3DS...
Wewnątrz lustrzanego 12-ścianu foremnego...

Wewnątrz torrusa...
Scena 3DS

Powierzchnia NURBS...

Przykładowa sceny składająca się z wyłącznie moich obiektów...



Powierzchnie NURBS....



Obiekt 3DS, około 60.000 trójkątów















Efekt rozszczepienia światła...




Powierzchnia NURBS – za mało trójkątów

Przykład rekurencyjnych odbić obiektu...


Efekt rozproszenia światła zastosowany na tle sceny...
Płytka szklana...



Kula...
Suszarka...


Obiekty wczytane z plików MD2 (QuakeII)

Przezroczysty torrus...
Obiekt 3DS (świątynia) około 60000 trójkątów...


		Plik 3DS bez interpolacji normalnych (trójkąty mają jedna normalna, obliczaną z iloczynu wektorowego boków)





















A oto wyrenderowane obrazy tej samej sceny dla trójkątów GPT (interpolacja w wierzchołkach, z różnymi właściwościami materiałowymi)









Wnętrze obiektu 3DS, 60000 trójkątów

























	1.	Bibliografia
[1] http://pl.wikipedia.org/wiki/Ray_Tracing   wikipedia, opis raytracingu
[2] http://gpurt.sourceforge.net/DA07_0405_Ray_Tracing_on_GPU-1.0.5.pdf GPU RT
[3]  – Forward Raytracing
[4] http://portal.acm.org/citation.cfm?id=1090144&dl=GUIDE&coll=GUIDE&CFID=68216571&CFTOKEN=23648247   Robust and numerically stable Bézier clipping method for ray tracing NURBS surfaces
[5] http://www.ocf.berkeley.edu/~jmich/res/cs184report.pdf
[6] http://www.cs.utah.edu/vissim/papers/raynurbs/node3.html – opis raytracingu NURBS
[7] http://pl.wikipedia.org/wiki/Forward_raytracing
[8] http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html – baza Bspline (używana przez powierzchnie NURBS)
[9] http://www.cs.utah.edu/vissim/papers/raynurbs/node4.html – podział płatków NURBS
[10] http://www.cs.utah.edu/vissim/papers/raynurbs/node6.html – rootfinder Newton,artefakty
[11]http://pl.wikipedia.org/wiki/Krzywa_B-sklejana
[12]https://www.cs.tcd.ie/courses/baict/bass/4ict10/Hillary2003/pdf/Lecture17_6Mar.pdf
porównanie metod lokalizacji: BVH, BSP, OctTree, Cells




















	Oświadczam, że pracę dyplomoiwą wykonałem samodzielnie.
	Łukasz Gryglicki

