# DeepTernary ile B ve C problarının karşılaştırılması

21 Kasım 2026. Bu protokol model sonuçları oluşmadan yazıldı. Amaç aynı iki bağlanma ucunu bilen bağımsız bir programın, iki protein arasındaki göreli yönelimi tahmin ettiğinde, tek T-box'ı örten çözümleri ilk sıralara getirip getirmediğini ölçmek.

## Tahmin

- Her cep için kimyasal bağları ve stereokimyası korunmuş tam prob kullanılır. B probunda kol 5QS1 fragmanının C4 atomundan, C probunda 5QSE fragmanının C6 atomundan çıkar. Cebin konumu ve CypA'ya bağlanan ucun konumu deneysel ikili yapılardan gelir. Geçmiş üçlü poz, DNA koordinatları, MD yörüngeleri ve örtüşme etiketleri tahmine verilmez.
- İki proteinin başlangıçtaki göreli pozunu taşımamak için bağımsız koordinat çerçeveleri kullanılır. Tam ligandın üç boyutlu konformasyonu bağ grafiğinden üretilir. Protein zincirleri ile ligand adları/atom eşlemeleri, giriş dosyaları, kaynak kimlikleri ve varsa çalışma zamanı uyumluluk değişiklikleri saklanır.
- DeepTernary'nin iki uçlu ligand modu, CPU'da aynı ağırlıklar, seçenekler ve her prob için aynı 40 rastgele tohumla kullanılır. 5T35/MZ1 resmî örneği kurulumu ve çıktıların okunmasını sınar. Bilinen bir örnektir ve TBXT üzerinde genelleme kanıtı sağlamaz. Modelin kendi tahmini protein 2 RMSD puanı küçükten büyüğe sıralanır. Tam eşitlikte küçük tohum önce gelir. DNA ölçümleri sıralamayı etkilemez.

## Sonradan ölçüm

- Başarı için aynı modelde, iki ucun da kendi proteini hizalandığında deneysel ikili yerleşimine en fazla 3 Å ağır atom RMSD ile oturması, kabul edilebilir sterik geometri ve tek T-box ile CypA örtüşmesi aranır. 2 Å ve 5 Å sınırları duyarlılık kontrolüdür. Proteinler arasındaki ya da probun herhangi bir proteinle yaptığı fiziksel çakışma ile referans DNA'yla aranan örtüşme ayrı kaydedilir. Kimyasal eşdeğer atomlar, çözülmemiş eşleme ve ligand geometrisi sapmaları ayrıca belirtilir.
- Her modelde TBXT'ye göre 8CDN hizalanır. Önceki çalışmayla aynı ölçüt kullanılır: uç nükleotitler çıkarıldıktan sonra, en az iki farklı nükleotitten toplam en az üç DNA ağır atomu CypA'ya 2,5 Å'dan daha yakınsa tek T-box alanı örtülür. DNA simülasyon/tahmin girdisi değildir.
- Birinci modeli, ilk 5, ilk 10 ve tüm 40 modeldeki sayıları, iki ucun doğru oturduğu modelleri ve varsa ilk ortak başarının sırasını raporlarız. Proteinler arası göreli yönelim deneysel olarak bilinmediğinden TBXT için model puanı deneysel RMSD veya afinite olarak yorumlanamaz. Aynı modele ait birden fazla tohum bağımsız bir biyolojik deney değildir.

## Çıktılar görülmeden sabitlenen uygulama ayrıntıları

İki protein arasındaki en kısa ağır atom uzaklığı en az 2,6 Å, her protein ile tam prob arasındaki en kısa ağır atom uzaklığı en az 2,2 Å aranır. B yapı girdisi 5QS1'den, C girdisi 5QSE'den gelir. C tarafında deneysel fragman pozunu koruyarak Asp177 yan zinciri çıkarılıp yabanıl tip Gly177 temsil edilir. Dolayısıyla sonuç, her bir prob ve ona eşlik eden yerel cep konformasyonu içindir.

DeepTernary'nin örnek komutu iki bağlanma ucunun atomlarını dosyalardaki uzaysal yakınlıktan tanımaya çalışır. Eski üçlü koordinatları kullanmadan çalışabilmek için B ve C problarında atom eşlemesi bağ grafiğinden sabitlenir. Her ikili proteinin cebindeki deneysel koordinatlar kendi çerçevesinde kalır. Liganda ait koordinatlar bağımsız olarak üretilen konformerin çerçevesindedir. Kullanılan giriş uyarlaması resmî kodun çalışır sürümü ve kontrol sonuçlarıyla birlikte arşivlenir.
