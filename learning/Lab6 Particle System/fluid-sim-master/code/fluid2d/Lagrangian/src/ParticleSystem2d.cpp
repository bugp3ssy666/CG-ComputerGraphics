#include "ParticleSystem2d.h"
#include "Lagrangian/include/ParticleSystem2d.h"
#include <iostream>
#include "Global.h"
#include <unordered_set>
#include <algorithm>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        ParticleSystem2d::ParticleSystem2d()
        {
        }

        ParticleSystem2d::~ParticleSystem2d()
        {
        }


        void ParticleSystem2d::setContainerSize(glm::vec2 lower = glm::vec2(-1.0f, -1.0f), glm::vec2 upper = glm::vec2(1.0f, 1.0f))
        {
            lower *=  Lagrangian2dPara::scale;
            upper *= Lagrangian2dPara::scale;

            lowerBound = lower - supportRadius + particleDiameter;
            upperBound = upper + supportRadius - particleDiameter;
            containerCenter = (lowerBound + upperBound) / 2.0f;

            glm::vec2 size = upperBound - lowerBound;

            blockNum.x = floor(size.x / supportRadius);
            blockNum.y = floor(size.y / supportRadius);

            blockSize = glm::vec2(size.x / blockNum.x, size.y / blockNum.y);

            blockIdOffs.resize(9);
            int p = 0;
            for (int j = -1; j <= 1; j++)
            {
                for (int i = -1; i <= 1; i++)
                {
                    blockIdOffs[p] = blockNum.x * j + i;
                    p++;
                }
            }

            particles.clear();
        }


        int ParticleSystem2d::addFluidBlock(glm::vec2 lowerCorner, glm::vec2 upperCorner, glm::vec2 v0, float particleSpace)
        {

            lowerCorner *= Lagrangian2dPara::scale;
            upperCorner *= Lagrangian2dPara::scale;

            glm::vec2 size = upperCorner - lowerCorner;

            if (lowerCorner.x < lowerBound.x ||
                lowerCorner.y < lowerBound.y ||
                upperCorner.x > upperBound.x ||
                upperCorner.y > upperBound.y)
            {
                return 0;
            }

            glm::uvec2 particleNum = glm::uvec2(size.x / particleSpace, size.y / particleSpace);
            std::vector<ParticleInfo2d> tempParticles(particleNum.x * particleNum.y);

            Glb::RandomGenerator rand;
            int p = 0;
            for (int idX = 0; idX < particleNum.x; idX++)
            {
                for (int idY = 0; idY < particleNum.y; idY++)
                {
                    float x = (idX + rand.GetUniformRandom()) * particleSpace;
                    float y = (idY + rand.GetUniformRandom()) * particleSpace;

                    tempParticles[p].position = lowerCorner + glm::vec2(x, y);
                    tempParticles[p].blockId = getBlockIdByPosition(tempParticles[p].position);
                    tempParticles[p].density = Lagrangian2dPara::density;
                    tempParticles[p].velocity = v0;
                    p++;
                }
            }

            particles.insert(particles.end(), tempParticles.begin(), tempParticles.end());
            return particles.size();
        }


        uint32_t ParticleSystem2d::getBlockIdByPosition(glm::vec2 position)
        {
            if (position.x < lowerBound.x ||
                position.y < lowerBound.y ||
                position.x > upperBound.x ||
                position.y > upperBound.y)
            {
                return -1;
            }

            glm::vec2 deltePos = position - lowerBound;
            uint32_t c = floor(deltePos.x / blockSize.x);
            uint32_t r = floor(deltePos.y / blockSize.y);
            return r * blockNum.x + c;
        }


        void ParticleSystem2d::updateBlockInfo()
        {
            std::sort(particles.begin(), particles.end(),
                      [=](ParticleInfo2d &first, ParticleInfo2d &second)
                      {
                          return first.blockId < second.blockId;
                      });

            blockExtens = std::vector<glm::uvec2>(blockNum.x * blockNum.y, glm::uvec2(0, 0));
            int curBlockId = 0;
            int left = 0;
            int right;
            for (right = 0; right < particles.size(); right++)
            {
                if (particles[right].blockId != curBlockId)
                {
                    blockExtens[curBlockId] = glm::uvec2(left, right);
                    left = right;
                    curBlockId = particles[right].blockId;
                }
            }
            blockExtens[curBlockId] = glm::uvec2(left, right);
        }

        // ========== 喷泉粒子生成 ==========
        // 在喷泉源位置生成新粒子，模拟向上喷出的水流
        // sourcePos: 喷泉源位置（归一化坐标）
        // velocity: 喷出速度
        // numParticles: 本次生成的粒子数量
        // spread: 水平扩散范围
        void ParticleSystem2d::addFountainParticles(glm::vec2 sourcePos, glm::vec2 velocity, int numParticles, float spread)
        {
            // 检查是否达到最大粒子数
            if (particles.size() >= maxParticles) return;

            // 将归一化坐标转换为世界坐标
            glm::vec2 worldPos = sourcePos * scale;
            
            Glb::RandomGenerator rand;

            for (int i = 0; i < numParticles; i++)
            {
                // 达到最大粒子数时停止
                if (particles.size() >= maxParticles) break;

                ParticleInfo2d newParticle;
                
                // 在源位置附近随机生成粒子位置
                float offsetX = (rand.GetUniformRandom() - 0.5f) * spread * scale;
                float offsetY = rand.GetUniformRandom() * particleDiameter;
                newParticle.position = worldPos + glm::vec2(offsetX, offsetY);

                // 检查位置是否在容器内
                if (newParticle.position.x < lowerBound.x ||
                    newParticle.position.x > upperBound.x ||
                    newParticle.position.y < lowerBound.y ||
                    newParticle.position.y > upperBound.y)
                {
                    continue;
                }

                // 设置初始速度，带有轻微随机扰动
                float velRandX = (rand.GetUniformRandom() - 0.5f) * 0.2f;
                float velRandY = rand.GetUniformRandom() * 0.1f;
                newParticle.velocity = velocity + glm::vec2(velRandX, velRandY);

                // 初始化其他属性
                newParticle.accleration = glm::vec2(0.0f);
                newParticle.density = Lagrangian2dPara::density;
                newParticle.pressure = 0.0f;
                newParticle.pressDivDens2 = 0.0f;
                newParticle.blockId = getBlockIdByPosition(newParticle.position);

                particles.push_back(newParticle);
            }
        }

        // ========== 移除死亡粒子 ==========
        // 移除超出容器边界的粒子，回收资源
        void ParticleSystem2d::removeDeadParticles()
        {
            // 使用erase-remove idiom移除无效粒子
            particles.erase(
                std::remove_if(particles.begin(), particles.end(),
                    [this](const ParticleInfo2d& p) {
                        // 如果粒子位置超出边界太多，标记为死亡
                        float margin = supportRadius * 2.0f;
                        return p.position.x < lowerBound.x - margin ||
                               p.position.x > upperBound.x + margin ||
                               p.position.y < lowerBound.y - margin ||
                               p.position.y > upperBound.y + margin ||
                               p.blockId == (uint32_t)-1;
                    }),
                particles.end()
            );
        }

    }
}
