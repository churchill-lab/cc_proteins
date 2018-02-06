docker network create --attachable qtlviewertest

docker run --rm -p 8000:8000 --name SOMENAME --network qtlviewertest -v /Users/mvincent/work/qtlviewer/data/Attie_DO378_eQTL_viewer_v3.RData:/app/qtlapi/data/ppppppppppppppdata.RData -v /Users/mvincent/work/qtlviewer/data/ccfoundersnps.sqlite:/app/qtlapi/data/ccfounders.sqlite -v /Users/mvincent/work/qtlviewer/qtlapi/qtlapi.R:/app/qtlapi/qtlapi.R mattjvincent/qtlapi /app/qtlapi/qtlapi.R

Than in docker-compose.yml...

networks:
  default:
    external:
      name: qtlviewertest
