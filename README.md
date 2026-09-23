# Центральность в BMN DST (run8 25.09)

Скрипт вычисляет центральность событий по DST-файлам BMN, применяет поправки на вершину и RefMult, отбор физических ранов, pile-up и триггер, и сохраняет результат в отдельный ROOT-файл. Полученный файл можно подключать как "friend chain" в своих анализах.

Главный скрипт `test_cent_bmndst.C`

## Что делает

Для каждого события вычисляет:

- **vertex-коррекции** (`vtxXcorr`, `vtxYcorr`, `vtxZcorr`) по run-by-run калибровкам;
- **RefMult** - число треков в кинематическом окне: `0.05 < pT < 2.0 GeV`, `0.7 < η < 2.5`, `DCA_R < 1 см`;
- **RefMultCorr** - RefMult с run-dependent поправкой;
- **centrality_Bin** - центральность (0…5 для 0–60%, `-1` иначе);
- **passSel** - флаг прохождения отбора события;
- **noPileup** - флаг отсутствия pile-up (метод O. Golosov, S. Sedyk).

## Требования

- Файлы поправок:
  - `run8_25.09_corr_VtxXYZ.root` - калибровки вершины;
  - `run8_25.09_corrections_DifRefRunId.root` - поправки RefMult по runId.

## Запуск

```cpp
// в ROOT
root -l -q -b test_cent_bmndst.C(
    "/path/to/mpd_run_Top_8130_ev0_p16.root",   // входной DST
    "./",                                        // директория выходных файлов
    "./run8_25.09_corr_VtxXYZ.root",            // калибровки вершины
    "./run8_25.09_corrections_DifRefRunId.root" // поправки RefMult
)

На выходе - файл вида: "cent_" + Имя_DST_файлам. Пример имени выходного файла для DST-файла `mpd_run_Top_8130_ev0_p16.root`: `cent_mpd_run_Top_8130_ev0_p16.root`

### файл `cent_mpd_run_Top_8130_ev0_p16.root` содержит дерево `bmndata` с ветками

| Ветка | Тип | Описание |
|---|---|---|
| `runId` | `UInt_t` | Номер рана |
| `evtId` | — | ID события |
| `noPileup` | `int` | `1` — нет pile-up, `0` — есть, `-1` — не определено |
| `PhysRunCut` | `bool` | Прошёл ли ран в список физических |
| `CCT2_mask` | — | Триггер CCT2 |
| `passSel` | `bool` | `1` — событие прошло все отборы |
| `RefMult` | `unsigned long` | Множественность |
| `RefMultCorr` | `float` | Скорректированная множественность |
| `centrality_Bin` | `int` | Центральность |

### Значения `centrality_Bin`

| Значение | Диапазон |
|---|---|
| `0` | 0–10% |
| `1` | 10–20% |
| `2` | 20–30% |
| `3` | 30–40% |
| `4` | 40–50% |
| `5` | 50–60% |
| `-1` | >60%, или событие не прошло отбор |

## Как использовать:

Подключите полученный файл как friend к исходному DST:

```cpp
TChain *chain = new TChain("bmndata");
chain->Add("mpd_run_Top_8130_ev0_p16.root");

TChain *chainCent = new TChain("bmndata");
chainCent->Add("cent_mpd_run_Top_8130_ev0_p16.root");

chain->AddFriend(chainCent);

Пример использования с ROOT::RDataFrame (см. функцию `UsageExample()`)

ROOT::RDataFrame dd(*chain);
auto h = dd.Filter("passSel==1").Histo1D({"h", ";cent", 6, -0.5, 5.5}, "centrality_Bin");
