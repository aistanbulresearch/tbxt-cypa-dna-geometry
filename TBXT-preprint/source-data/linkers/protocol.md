# B cebi bağlantı kolu karşılaştırması

21 Kasım 2026. Bir TBXT zinciri (6F58 A), aynı modellenmiş 112–119 döngüsü, 5QS1 NZ1/C4 ve 9GHY A1ILH/NBR bağlanma uçları. Ar, B cebindeki aromatik karbonu; R, CypA'ya tutunan ucu gösterir. Dört nötr trans-amid prob:

| Kod | Ar ve NH–R arasındaki kol |
|---|---|
| ref | CH2–O–CH2–CH2–O–CH2–C(=O) |
| short | CH2–O–CH2–C(=O) |
| long | CH2–O–CH2–CH2–O–CH2–CH2–O–CH2–C(=O) |
| trans | CH2–O–CH2–CH=CH–CH2–C(=O), E izomeri |

Son analog referansla aynı sayıda omurga atomuna sahiptir; hem esnekliği hem oksijen sayısı değişir. Kısa ve uzun analoglarda bir oksietilen birimi çıkarılır veya eklenir.

Önceden tanımlanmış geometrik tarama: her kol için aynı dört tohum (98101–98104) ve tohum başına 18.000 DNA-kör burulma örneği. Kabul için bağ dışı çakışma filtreleri, CypA–TBXT ağır atomları arasında en az 2,6 Å; 3,0 Å duyarlılık kontrolü. DNA'yı yalnızca bu filtrelerden sonra 8CDN tek T-box referansı olarak değerlendir. İki uç nükleotidi her zincirde dışarıda bırakmak için ilk ve son nükleotidi çıkar. 2,5 Å içinde en az üç ağır atom ve en az iki farklı nükleotit örtüşme sayılır. Sayılar bir arama algoritmasının sonuçlarıdır; fiziksel bağlanma olasılığı değildir.

Hareket karşılaştırmasında her kimya için ilk örtüşmeyen yakın (4–8 Å) ve uzak (10–18 Å) başlangıç ve varsa ilk örten başlangıç kullanılacak. İlk giriş, kalış, son kare ve iki bağlanma ucunun temasları ayrı kaydedilecek. Başlangıç aralıkları tek T-box iç DNA'sına mesafedir. İki ucun başlangıçta bağlı olduğu sistemlerin sonuçları çözeltide serbest moleküllerin bir araya gelme olasılığını vermez. Koşullar veya örnekleme yetersizse üstün kol ilan edilmeyecek.

Geometrik taramadan sonra, dinamik sonuçlara bakılmadan belirlenen uygulama: çakışmasız yakın, uzak ve örten birer başlangıcı olan mevcut, uzun ve E-alken kollarında iki hız tohumu (33512, 33513), toplam 18 çözelti MD koşusu. Ana çalışmayla aynı TIP3P su, 0,15 M tuz, 300 K, 1 bar, ff14SB protein ve Sage 2.2.1/MMFF94 yüklü prob, kısıtlı 4 ps dengeleme ve kısıtsız 50 ps üretim. Birinci yakın başlangıç serisi tüm kimyalarda bir CPU iş parçacığında yürütülebilir. Kalan koşular GPU'da yürütülür. Kullanılan cihaz her sonuçta kaydedilir. Koşular arası karşılaştırmalar aynı başlangıç türü ve hız tohumu içinde yapılır. Bir erken geometrik kontrol için iki tekrar kullanıyoruz; bunlar denge olasılığını tahmin etmek için yeterli olmayabilir.

Örtüşmeye ek olarak, her iki prob ucunun kendi başlangıç protein temaslarının en az %50'sini hâlâ koruduğu örten kareleri say. Başlangıç ve karelerdeki temas mesafesi 4,5 Å.

İlk taramada geçerli yerleşimlerin seyrek olması ve kısa kolda hiç çıkmaması nedeniyle, ayrı bir örnekleme kontrolü eklendi. Bu karar ilk 16 tohum grubunun sonuçlarından sonra, ek tohumlar hesaplanmadan önce alındı. Dört kolun her biri için aynı 24 yeni tohum (99201–99224) × 18.000 deneme yapılır. Aynı ağır atom filtreleri ve iki ayrı payda (tüm denemeler; sığan konumlar) korunur. Sonuçlar [ek tarama tablosunda](holdout-screen.csv) verilir. Ana MD'nin önceden seçilen başlangıçları korunur. Yeni kısa kol konumları bulunursa ayrıca raporlanır.

24 ek tohumun ardından bulunan tek kısa kol konumu, başlangıçta 9,64 Å'dadır. İki mevcut hız tohumu ile bu ayrı konumda iki ilave 50 ps MD koşusu planlandı. Yakınlık kategorileri dışına düştüğü için kısa kol, üç kollu eşleştirilmiş 18 tekrarın sıralamasına eklenmez. Ondan gelen sonuçlar bağımsız bir duyarlılık kontrolü olarak değerlendirilir.
