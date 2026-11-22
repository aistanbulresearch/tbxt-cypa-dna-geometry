# Stereokimya kontrolü ve yeniden hesaplama kararı

22 Kasım 2026. İlk 40 B ve 40 C tahmini tamamlandıktan sonra stereokimya denetiminde B probunun iki karbon merkezinin (CAJ ve CBK) ters olduğu görüldü. Aynı örnekleme birkaç tohumla yeniden çalıştırıldığında değişikliğin model tahmininden önce, tam prob için yeni bir konformer üretilirken oluştuğu doğrulandı. Başlangıç SDF'sinde tanımlanan stereoizomer doğruydu. İlk serinin puanları, modelleri ve ölçümleri karşılaştırılabilirlik için saklandı.

İstenen stereoizomerlerin karşılaştırılması için B ve C, ayrıca deneysel yapısı bilinen 5T35, aynı yöntemle yeniden çalıştırılır. ETKDGv3 örneklemesinde hidrojenler açıkça temsil edilir. Her deneysel uç için yalnızca kendi içindeki uzaklıklar kullanılır. İki ucun sınırları tek kez birlikte tutarlı hale getirilir. Her iki probda da ±0,10 Å sınırları bu tutarlılık denetiminden geçmediği için, ön denemede her ikisinin de geçtiği ±0,20 Å kullanılır. Uçların birbirine göre konumu serbesttir.

Her tohumda ağ değerlendirmesinden önce gerçek üç boyutlu koordinatlardan stereokimya yeniden tayin edilir. Önceden tanımlanmış kiral merkezlerin tamamı ve tanımlanmış çift bağ geometrileri girişteki kimyasal grafikle eşleşmelidir. Başarısız konformer üretimi ve deterministik tekrar tohumları kaydedilir. Ağın ürettiği ham ligand ve torsiyon düzeltmesinden sonraki ligand ayrıca denetlenir. Puan, DNA örtüşmesi ve iki proteinin tahmin edilen göreli pozu konformer seçiminde kullanılmaz.

Birincil bilimsel karşılaştırma bu denetimden geçen tam problarla yapılan 40'ar tahmindir. Sıralama, iki uç için konum sapması, sterik kabul ölçütleri ve sonradan yapılan tek T-box örtüşme ölçümü [önceden yazılan protokoldeki](Protokol.md) gibidir. İlk serinin ikinci koordinat çerçevesi deneyi tarihsel bir duyarlılık kontrolüdür.

İlk stereokimya denemesinde 5T35'in ilk tohumu 24 konformer üretim denemesini tüketti. Bunun üzerine üç sistemin tamamı için ETKDG'nin rastgele başlangıç koordinatları seçeneğine geçildi. Yarım kalan deneme ayrı saklanır. Bu kararın gerekçesi konformer üretiminin başarısıdır. Bu aşamada 5T35 deneysel üçlü pozu örnekleme veya tohum seçiminde kullanılmadı.
