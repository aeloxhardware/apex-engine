# Convention du Dev Tracker

_Extrait de l'ancien CLAUDE.md._

## Dev Tracker (convention de rédaction du suivi de développement)

Le « Dev Tracker » (« DT ») est la liste de suivi de ce qui a été développé et de
ce qui reste à faire. **À rédiger dans deux cas seulement** : quand c'est
nécessaire (fin d'un chantier, changement d'état important, reprise après une
interruption, décision de l'utilisateur qui déplace une priorité), ou quand
l'utilisateur le demande (« écris le DT »). Pas à chaque message.

**Marqueurs, un par ligne, toujours en tête :**
- 🟢 **résolu** : fait, ou question tranchée par l'utilisateur ;
- ⏳ **en cours** : tâche commencée et non terminée (au plus une ou deux à la fois) ;
- ⬛ **non traité, sans successeur** : rien ne dépend de ce point, il ne bloque rien ;
- ♦️ **à résoudre avant de continuer** : bloquant, c'est ce qui conditionne la suite.

**Structure (dans cet ordre) :**
1. *Ce qui a été développé* : 🟢 (correctifs, puis nouvelles fonctions) ;
2. *Réponses de l'utilisateur* : 🟢 ce qui est tranché ;
3. *À résoudre avant de continuer* : ♦️ ;
4. *Sans bloquer* : ⬛ ;
5. *Fait* : rappel condensé des 🟢 inchangés.

**Où il vit** : `DT.md` à la racine du projet (version courante, à mettre à jour à chaque
rédaction) ; `LAST_SESSION.md` en garde le contexte de session.

**Règles de rédaction :**
- **Vérifier avant de classer.** Un point n'est 🟢 que si le code ou une mesure le
  confirme, jamais sur la foi d'un document : `LAST_SESSION.md` a déjà annoncé
  « non commités » des correctifs qui l'étaient (vérifié le 2026-10-06).
- **Un point ♦️ dit pourquoi il bloque** et ce qui le débloquera.
- **Actualiser la liste existante** (changer les marqueurs, ajouter les nouveautés)
  plutôt que la réécrire d'un autre style.
- Citer les chiffres mesurés et les limites (#NN) plutôt que des généralités.

