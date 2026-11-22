# Tek brachyury zinciriyle DNA yüzeyinin örtülmesi

## Biyolojik soru

Brachyury bir transkripsiyon faktörüdür. Belirli DNA dizilerine bağlanıp genlerin ne ölçüde çalışacağını düzenler. Kordoma hücreleri bu programa bağımlıdır. İki işlevli bir probun bir ucunun brachyury'ye, diğer ucunun hücrenin siklofilin A (CypA) proteinine bağlandığını varsayarak CypA'nın DNA'nın gerektirdiği yeri kaplayıp kaplayamayacağını sınadık. Fragmanların B ve C ceplerine bağlanması kristalografide gözlenmiştir. Probun ve CypA'nın birlikte yerleşimi ise hesaplamalı bir modeldir. [Kordomada işlevsel genetik](https://www.nature.com/articles/s41591-018-0312-3) · [Brachyury kristal yapıları ve cepler](https://www.nature.com/articles/s41467-025-56213-1)

## Başlangıçların seçimi

Hedef, 6F58'deki yabanıl tip brachyury'nin A zinciridir. İkinci brachyury kopyası ve DNA hareketli sistemden çıkarılır. B'de 5QS1 fragmanının C4 atomundan, C'de 5QSE fragmanının C6 atomundan çıkan probun kimyası her cepte sabit tutulur. B ve C probunun brachyury'ye bağlanan uçları farklıdır.

5QS1, yabanıl tip brachyury'nin fragmanla kristal yapısıdır. 5QSE'nin fragman pozu ise G177D varyantında gözlenmiştir. Bu pozu yabanıl tip proteine aktarmak C adayı için ayrıca bir modelleme varsayımıdır. [5QS1](https://www.rcsb.org/structure/5QS1) · [5QSE](https://www.rcsb.org/structure/5QSE)

| Başlangıç | Anlamı |
|---|---|
| B-on / C-on | Önceki iki zincirli hesapta kullanılan ve tam 6F58 DNA'sını örten poz. Bunlardan yalnız B-on tek T-box DNA'sını da örter. |
| B-off / C-off | Önceki sterik taramada DNA'dan uzak seçilmiş yönelim. 8CDN iç DNA'sına B'de 14,12 Å, C'de 15,08 Å uzaklık. |
| B-near / C-near | Aynı taramanın sırasına göre seçilmiş, DNA'ya daha yakın ve başlangıçta örtmeyen ikinci yönelim. Aynı uzaklıklar B'de 7,51 Å, C'de 7,82 Å. |

Örtmeyen dört başlangıcın her birinde CypA ile 6F58 veya 8CDN'nin uçlar dahil herhangi bir DNA ağır atomu arasında 2,5 Å'dan kısa mesafe yoktu. Başlangıçlar, iki brachyury kopyasının da bulunduğu eski yapı taramasından seçildi. Yeni minimizasyon, dengeleme ve üretim hesaplarının tamamı tek kopyayla yapılıyor. Başlangıçları seçerken üretim yörüngeleri kullanılmadı. Her iki cep için 112–119 döngüsünün koordinatları 5QSE'den gelir. C modelindeki Arg54 de 5QSE'de gözlenen yan zincir konumunu kullanır.

## Moleküler dinamik ve ölçüm

Her başlangıç geometrisine üç hız tohumu atanır. Küçük molekülün Sage 2.2.1 parametreleri ve pilotta kullanılan MMFF94 kısmi yükleri, proteinin ff14SB kuvvet alanı, açık TIP3P su ve 0,15 M NaCl önceki hesapla aynı tutulur. Hedef sıcaklık 300 K, basınç 1 bar, zaman adımı 2 fs'dir. Minimizasyondan sonra 2+2 ps kademeli konum kısıtlaması ve 50 ps kısıtsız üretim uygulanır. Üretim 0–50 ps dahil 1 ps aralıkla kaydedilir. Aynı yörüngenin 51 karesi birbirine bağımlıdır.

Birincil ölçümde her karede brachyury'nin ortak çekirdek Cα atomları hizalanır. 8CDN'deki iki DNA zincirinin uçları çıkarıldıktan sonra, en az iki farklı nükleotitte toplam en az üç DNA ağır atomunun CypA'ya 2,5 Å'dan yakın olması örtüşme sayılır. Kaynak poz, minimizasyon, kısıtlı dengeleme ve kısıtsız üretim ayrı değerlendirilir. Üretime zaten örtüşerek giren bir koşu, üretim sırasında kendiliğinden geçiş yaptığı şeklinde yorumlanmaz. Eşlenmiş 6F58 DNA bölgesi ve tam 6F58 DNA ayrı ölçülür. [Tek T-box DNA yapısı](https://www.rcsb.org/structure/8CDN) · [Palindromik DNA yapısı](https://www.rcsb.org/structure/6F58)

## Duyarlılık analizi

Prob uçlarının yer değiştirmesine duyarlılık ayrıca incelendi. Bu ek ölçüt yörüngeler görüldükten sonra tanımlandı. Tek T-box örtüşmesinin olduğu karelerde, brachyury ucunun brachyury'ye ve CypA ucunun CypA'ya göre kaymasının aynı anda 2 Å veya 2,5 Å ağır atom RMSD sınırında kalması arandı. Sonuçlar [iki ucun poz kontrolü tablosunda](old-monomer-heads.csv) verilir.

Kısa ve önceden bağlanmış komplekslerden başlatılan bu pilot, çözeltide kompleksin kendiliğinden oluşmasını, DNA'nın yerinden ayrılmasını veya kordoma tedavisindeki etkinliği ölçmez.
