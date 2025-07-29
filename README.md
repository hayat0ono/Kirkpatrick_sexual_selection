# Kirkpatrickの性選択シミュレーション

このプロジェクトは、Kirkpatrickの性選択理論に基づいた個体群シミュレーションを実装したものです。

## 概要

このシミュレーションでは、以下の要素を考慮した性選択の進化をモデル化しています：

- **T遺伝子**: オスの表現型を決定する遺伝子（T0, T1）
- **P遺伝子**: メスの好みを決定する遺伝子（P0, P1）
- **選択係数**: メスの好みの強さを表すパラメータ
- **生存率**: 遺伝子型による生存確率の違い

## プロジェクト構造

```
Kirkpatrick_sexual_selection/
├── README.md
├── requirements.txt
├── .gitignore
└── sexual_selection/
    ├── model.py
    ├── male.py
    ├── female.py
    ├── individual.py
    └── run/
        ├── runaway.py
        └── equilibrium_state.py
```

## 主要なファイル

### モデル関連
- `sexual_selection/model.py`: メインのシミュレーションモデル
- `sexual_selection/male.py`: オス個体のクラス
- `sexual_selection/female.py`: メス個体のクラス
- `sexual_selection/individual.py`: 基本個体クラス

### 実行スクリプト
- `sexual_selection/run/equilibrium_state.py`: 平衡状態のシミュレーション
- `sexual_selection/run/runaway.py`: ランナウェイ効果のシミュレーション

## パラメータ説明

- **N**: 総個体数
- **male_female_ratio**: オスメス比（0.0-1.0）
- **t_gene_ratio**: T1遺伝子の初期頻度（0.0-1.0）
- **p_gene_ratio**: P1遺伝子の初期頻度（0.0-1.0）
- **a0_coeff**: 選択係数a0（Kirkpatrick理論のパラメータ）
- **a1_coeff**: 選択係数a1（Kirkpatrick理論のパラメータ）
- **s_ratio**: 生存率の違い（0.0-1.0）
- **male_death_prob**: オスの基本死亡確率
- **female_death_prob**: メスの死亡確率
- **mutation_rate**: 突然変異率（0.0-1.0）
- **lifetime**: 個体の寿命
- **num_child**: 1回の交配での子の数


## 必要なライブラリ

- Python 3.7+
- NumPy
- Matplotlib

## インストール

```bash
cd Kirkpatrick_sexual_selection
pip install numpy matplotlib
```


## 参考文献

- Kirkpatrick, M. (1982). Sexual selection and the evolution of female choice. Evolution, 36(1), 1-12.
- 伊庭斉志 (2007). 複雑系のシミュレーション - Swarmによるマルチエージェント・システム -. コロナ社.
- [東京大学伊庭研究室 - Swarmソフトウェアサポートページ](http://www.iba.t.u-tokyo.ac.jp/software/Swarm_Software/support_swarm.html)
